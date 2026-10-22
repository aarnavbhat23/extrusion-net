"""Complete-mixer timing and saved-activation scaling against a dense reference."""
import argparse
import json
import platform
import statistics
import time
from pathlib import Path

import torch

from models import ExtrusionMixer


def dense_forward(model, x):
    batch, n, d = x.shape
    weights = model.contacts(x)
    values = model.value(x).view(batch, n, model.heads, d // model.heads).transpose(1, 2)
    msg = (weights @ values).transpose(1, 2).reshape(batch, n, d)
    mass = torch.log1p(weights.sum(-1)).transpose(1, 2)
    return model.out(torch.cat([msg, mass], dim=-1))


def measure(n, mode, backward):
    torch.manual_seed(11)
    model = ExtrusionMixer(32, 2, saturate=False)
    x = torch.randn(1, n, 32, requires_grad=backward)
    function = model if mode == 'scan' else lambda value: dense_forward(model, value)
    saved = []
    def record(tensor):
        saved.append(tensor.numel() * tensor.element_size())
        return tensor
    with torch.autograd.graph.saved_tensors_hooks(record, lambda tensor: tensor):
        with torch.set_grad_enabled(backward):
            y = function(x)
            if backward:
                y.square().mean().backward()
    times = []
    for _ in range(3):
        model.zero_grad(set_to_none=True)
        x.grad = None
        start = time.perf_counter()
        with torch.set_grad_enabled(backward):
            result = function(x)
            if backward:
                result.square().mean().backward()
        times.append(time.perf_counter() - start)
    return {'length':n, 'mode':mode, 'backward':backward,
            'median_seconds':statistics.median(times), 'seconds':times,
            'saved_tensor_logical_bytes':sum(saved),
            'note':'Saved tensor bytes count logical views repeatedly; this is not peak resident memory.'}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', default='results/mixer_scaling.json')
    args = parser.parse_args()
    torch.set_num_threads(1)
    rows = []
    for backward in (False, True):
        for n in (128, 512, 1024, 4096):
            for mode in ('scan','dense'):
                if mode == 'dense' and n > 1024:
                    continue
                row = measure(n, mode, backward)
                rows.append(row)
                print(json.dumps(row), flush=True)
    result = {'torch':torch.__version__, 'python':platform.python_version(),
              'machine':platform.machine(), 'threads':1,'rows':rows}
    Path(args.out).write_text(json.dumps(result, indent=2)+'\n')


if __name__ == '__main__':
    main()
