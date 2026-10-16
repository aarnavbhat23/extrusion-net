"""Extrusion Networks and the controls they are compared against.

Every model shares one backbone (embedding, local convolution, MLP). Only the
long-range mixer differs, so any difference in results comes from the mixer.
"""

import math

import torch
import torch.nn as nn
import torch.nn.functional as F

NEG = -1e9  # finite stand-in for log(0); -inf breaks logcumsumexp gradients
MAX_STOP_LOGIT = 12.0  # a barrier lets at most ~6e-6 of arms through


def _shift(x):
    """x[..., t] -> x[..., t-1], with 0 for t = 0."""
    return torch.cat([torch.zeros_like(x[..., :1]), x[..., :-1]], -1)


def extrusion_log_contacts(log_load, a, b):
    """log D(i, j): chance one motor holds i and j in the same loop.

    A motor loads on the bond just left of position k, then a left arm and a
    right arm move outward. a[m] = log P(left arm passes m) and
    b[m] = log P(right arm passes m). For i < j:

        D(i, j) = sum_{k=i+1..j} load_k * prod_{m=i+1..k-1} pass_left_m
                                        * prod_{m=k..j-1}   pass_right_m

    Computed exactly in O(n^2) with prefix sums. Shapes: (..., n) -> (..., n, n).
    The diagonal is NEG (self-contact is handled by the residual stream).
    """
    n = a.shape[-1]
    A, B = a.cumsum(-1), b.cumsum(-1)
    c = log_load + _shift(A) - _shift(B)  # term for each load bond k
    idx = torch.arange(n, device=a.device)
    after_i = idx[None, :] > idx[:, None]  # [i, k]: k >= i + 1
    C = c.unsqueeze(-2).expand(*c.shape[:-1], n, n).masked_fill(~after_i, NEG)
    R = torch.logcumsumexp(C, dim=-1)  # R[i, j] = logsumexp_{k=i+1..j} c_k
    upper = _shift(B).unsqueeze(-2) - A.unsqueeze(-1) + R
    upper = upper.masked_fill(~after_i, NEG)
    return torch.maximum(upper, upper.transpose(-1, -2))


def insulation_log_contacts(a, b):
    """Ablation: one arm launched from the query itself, stopped by barriers.

    Right of i: prod_{m=i+1..j-1} pass_right_m. Left of i: prod_{m=j+1..i-1}
    pass_left_m. This is close to stick-breaking attention with per-position
    stops; it lacks the two-arm, load-anywhere mechanism.
    """
    n = a.shape[-1]
    A, B = a.cumsum(-1), b.cumsum(-1)
    idx = torch.arange(n, device=a.device)
    right = _shift(B).unsqueeze(-2) - B.unsqueeze(-1)  # [i, j] = B_{j-1} - B_i
    left = _shift(A).unsqueeze(-1) - A.unsqueeze(-2)  # [i, j] = A_{i-1} - A_j
    j_gt_i = idx[None, :] > idx[:, None]
    j_lt_i = idx[None, :] < idx[:, None]
    return torch.where(j_gt_i, right, torch.where(j_lt_i, left, torch.full_like(right, NEG)))


class ExtrusionMixer(nn.Module):
    """Positions talk only through contacts made by simulated loop extrusion."""

    def __init__(self, d, heads, two_arm=True):
        super().__init__()
        self.heads, self.two_arm = heads, two_arm
        self.control = nn.Linear(d, 3 * heads)  # left stop, right stop, load
        self.value = nn.Linear(d, d)
        self.out = nn.Linear(d + heads, d)
        self.beta_raw = nn.Parameter(torch.zeros(heads))
        with torch.no_grad():
            self.control.bias.view(heads, 3)[:, :2] = -2.0  # leaky barriers at start

    def contacts(self, x):
        bsz, n, _ = x.shape
        z = self.control(x).view(bsz, n, self.heads, 3).permute(0, 2, 3, 1)
        a = -F.softplus(z[:, :, 0].clamp(max=MAX_STOP_LOGIT))  # log(1 - s_left)
        b = -F.softplus(z[:, :, 1].clamp(max=MAX_STOP_LOGIT))  # log(1 - s_right)
        if self.two_arm:
            log_d = extrusion_log_contacts(F.logsigmoid(z[:, :, 2]), a, b)
        else:
            log_d = insulation_log_contacts(a, b)
        beta = F.softplus(self.beta_raw)[None, :, None, None] * 4.0
        return -torch.expm1(-beta * log_d.exp())  # Poisson: P(at least one motor)

    def forward(self, x):
        bsz, n, d = x.shape
        w = self.contacts(x)
        v = self.value(x).view(bsz, n, self.heads, d // self.heads).transpose(1, 2)
        msg = (w @ v).transpose(1, 2).reshape(bsz, n, d)
        mass = torch.log1p(w.sum(-1)).transpose(1, 2)
        return self.out(torch.cat([msg, mass], -1))


class AttentionMixer(nn.Module):
    """Bidirectional softmax attention with RoPE or ALiBi positions."""

    def __init__(self, d, heads, pos):
        super().__init__()
        self.heads, self.pos = heads, pos
        self.qkv = nn.Linear(d, 3 * d)
        self.out = nn.Linear(d, d)
        slopes = torch.tensor([2 ** (-8 * (h + 1) / heads) for h in range(heads)])
        self.register_buffer("slopes", slopes)

    def _rope(self, x):
        n, dh = x.shape[-2], x.shape[-1]
        freqs = 1.0 / (10000 ** (torch.arange(0, dh, 2, device=x.device) / dh))
        ang = torch.arange(n, device=x.device)[:, None] * freqs[None, :]
        cos, sin = ang.cos(), ang.sin()
        x1, x2 = x[..., 0::2], x[..., 1::2]
        return torch.stack([x1 * cos - x2 * sin, x1 * sin + x2 * cos], -1).flatten(-2)

    def forward(self, x):
        bsz, n, d = x.shape
        dh = d // self.heads
        q, k, v = self.qkv(x).view(bsz, n, 3, self.heads, dh).permute(2, 0, 3, 1, 4)
        bias = None
        if self.pos == "rope":
            q, k = self._rope(q), self._rope(k)
        else:
            idx = torch.arange(n, device=x.device)
            dist = (idx[None, :] - idx[:, None]).abs().float()
            bias = -self.slopes[:, None, None] * dist
        y = F.scaled_dot_product_attention(q, k, v, attn_mask=bias)
        return self.out(y.transpose(1, 2).reshape(bsz, n, d))


class LSTMMixer(nn.Module):
    def __init__(self, d):
        super().__init__()
        self.lstm = nn.LSTM(d, d // 2, batch_first=True, bidirectional=True)

    def forward(self, x):
        return self.lstm(x)[0]


class Block(nn.Module):
    def __init__(self, d, mixer):
        super().__init__()
        self.ln_conv, self.ln_mix, self.ln_mlp = nn.LayerNorm(d), nn.LayerNorm(d), nn.LayerNorm(d)
        self.conv = nn.Conv1d(d, d, 5, padding=2)
        self.mixer = mixer
        self.mlp = nn.Sequential(nn.Linear(d, 2 * d), nn.GELU(), nn.Linear(2 * d, d))

    def forward(self, h):
        h = h + F.gelu(self.conv(self.ln_conv(h).transpose(1, 2)).transpose(1, 2))
        h = h + self.mixer(self.ln_mix(h))
        return h + self.mlp(self.ln_mlp(h))


MIXERS = {
    "extrusion": lambda d, H: ExtrusionMixer(d, H, two_arm=True),
    "one_arm": lambda d, H: ExtrusionMixer(d, H, two_arm=False),
    "attn_rope": lambda d, H: AttentionMixer(d, H, "rope"),
    "attn_alibi": lambda d, H: AttentionMixer(d, H, "alibi"),
    "bilstm": lambda d, H: LSTMMixer(d),
}


class SequenceModel(nn.Module):
    def __init__(self, mixer, vocab=7, n_classes=4, d=64, heads=4, layers=3):
        super().__init__()
        self.embed = nn.Embedding(vocab, d)
        self.blocks = nn.ModuleList(Block(d, MIXERS[mixer](d, heads)) for _ in range(layers))
        self.ln = nn.LayerNorm(d)
        self.head = nn.Linear(d, n_classes)

    def forward(self, tokens):
        h = self.embed(tokens)
        for block in self.blocks:
            h = block(h)
        return self.head(self.ln(h))


def n_params(model):
    return sum(p.numel() for p in model.parameters())
