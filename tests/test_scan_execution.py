"""Verify the scan's end-to-end mathematical and memory contracts."""
import copy
import unittest
from unittest.mock import patch

import numpy as np
import torch

from models import CascadedScanMixer, ExtrusionMixer, SequenceModel, extrusion_linear_scan
from test_models import brute_extrusion


def dense_forward(mixer, x):
    batch, n, width = x.shape
    weights = mixer.contacts(x)
    values = mixer.value(x).view(batch, n, mixer.heads, width // mixer.heads).transpose(1, 2)
    messages = (weights @ values).transpose(1, 2).reshape(batch, n, width)
    mass = torch.log1p(weights.sum(-1)).transpose(1, 2)
    return mixer.out(torch.cat([messages, mass], dim=-1))


class TestScanExecution(unittest.TestCase):
    def test_mixer_outputs_and_all_gradients_match_dense(self):
        torch.manual_seed(17)
        scan = ExtrusionMixer(12, 3, saturate=False).double()
        dense = copy.deepcopy(scan)
        x = torch.randn(2, 17, 12, dtype=torch.float64, requires_grad=True)
        reference_x = x.detach().clone().requires_grad_()
        actual = scan(x)
        expected = dense_forward(dense, reference_x)
        torch.testing.assert_close(actual, expected, rtol=1e-10, atol=1e-11)
        probe = torch.randn_like(actual)
        (actual * probe).sum().backward()
        (expected * probe).sum().backward()
        torch.testing.assert_close(x.grad, reference_x.grad, rtol=1e-9, atol=1e-10)
        for (name, parameter), (_, reference) in zip(scan.named_parameters(), dense.named_parameters()):
            self.assertIsNotNone(parameter.grad, name)
            torch.testing.assert_close(parameter.grad, reference.grad, rtol=1e-9, atol=1e-10, msg=name)

    def test_scan_handles_exact_barriers_and_singletons(self):
        rng = np.random.default_rng(7)
        for n in (1, 2, 13):
            load = rng.uniform(0, 1, n)
            pl, pr = rng.integers(0, 2, (2, n)).astype(float)
            values = rng.normal(size=(n, 4))
            result = extrusion_linear_scan(*(torch.tensor(a) for a in (load, pl, pr, values)))
            np.testing.assert_allclose(result.numpy(), brute_extrusion(load, pl, pr) @ values,
                                       rtol=1e-12, atol=1e-12)

    def test_untied_cascade_contains_linear_extrusion(self):
        torch.manual_seed(31)
        physical = ExtrusionMixer(12, 3, saturate=False).double()
        control = CascadedScanMixer(12, 3).double()
        state = physical.state_dict()
        for name in ('value.weight', 'value.bias', 'out.weight', 'out.bias', 'beta_raw'):
            control.state_dict()[name].copy_(state[name])
        with torch.no_grad():
            control.control.weight.copy_(physical.control.weight.view(3, 3, 12).repeat(1, 2, 1).reshape(18, 12))
            control.control.bias.copy_(physical.control.bias.view(3, 3).repeat(1, 2).reshape(18))
        inputs = torch.randn(2, 19, 12, dtype=torch.float64)
        torch.testing.assert_close(control(inputs), physical(inputs), rtol=1e-12, atol=1e-12)

    def test_full_linear_model_never_materializes_contacts(self):
        model = SequenceModel('extrusion_linear', d=16, heads=2, layers=2)
        with patch('models.extrusion_log_contacts', side_effect=AssertionError('dense contacts called')):
            model(torch.randint(0, 7, (1, 257))).square().mean().backward()

    def test_saved_activations_grow_linearly(self):
        totals = []
        for n in (128, 256):
            mixer = ExtrusionMixer(16, 2, saturate=False)
            sizes = []
            def record(tensor):
                self.assertNotEqual(tuple(tensor.shape[-2:]), (n, n))
                sizes.append(tensor.numel())
                return tensor
            with torch.autograd.graph.saved_tensors_hooks(record, lambda tensor: tensor):
                mixer(torch.randn(1, n, 16, requires_grad=True)).square().mean().backward()
            totals.append(sum(sizes))
        self.assertLess(totals[1] / totals[0], 2.2)


if __name__ == '__main__':
    unittest.main()
