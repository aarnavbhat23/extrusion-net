import unittest
import torch
from moment_text import MomentLanguageModel


class MomentTextTests(unittest.TestCase):
    def test_full_text_models_cannot_read_future_characters(self):
        torch.manual_seed(98)
        for kind in ('atom','mean','static','attention'):
            model=MomentLanguageModel(kind,32).double().eval();tokens=torch.randint(32,(2,16));changed=tokens.clone();changed[:,8:]=torch.randint(32,(2,8))
            first=model(tokens);second=model(changed)
            torch.testing.assert_close(first[:,:8],second[:,:8],rtol=1e-7,atol=1e-7)
            first.square().mean().backward()
            self.assertTrue(all(p.grad is None or torch.isfinite(p.grad).all() for p in model.parameters()))


if __name__=='__main__':unittest.main()
