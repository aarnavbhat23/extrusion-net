import unittest
import torch
from moment_memory import moment_state,atom_scores,value_codes,RetrievalOperator


class MomentMemoryTests(unittest.TestCase):
    def test_atom_condition_contains_all_written_values_and_matches_psd_test(self):
        torch.manual_seed(92);keys=torch.softmax(torch.randn(2,9,7,dtype=torch.float64)*2,dim=-1)
        keys=torch.nn.functional.normalize(keys,dim=-1);codes=value_codes(dtype=torch.float64)
        labels=torch.randint(8,(2,9));values=codes[labels];state=moment_state(keys,values)
        scores=atom_scores(state,keys,codes)
        self.assertTrue(torch.all(scores.gather(-1,labels[...,None])<=1+1e-10))
        gram=torch.einsum('bqd,bdij->bqij',keys,state)+1e-3*torch.eye(4,dtype=torch.float64)
        augmented=torch.cat([torch.ones(8,1,dtype=torch.float64),codes],dim=1)
        residual=gram[:,:,None]-augmented[None,None,:,:,None]*augmented[None,None,:,None,:]
        psd=torch.linalg.eigvalsh(residual)[...,0]>=-1e-10
        self.assertTrue(torch.equal(psd,scores<=1+1e-10))

    def test_inverse_moment_read_gradient_matches_finite_differences(self):
        torch.manual_seed(94)
        keys=(torch.rand(1,3,2,dtype=torch.float64)+.5).requires_grad_()
        queries=(torch.rand(1,2,2,dtype=torch.float64)+.5).requires_grad_()
        values=value_codes(dtype=torch.float64)[torch.tensor([[0,2,7]])]
        codes=value_codes(dtype=torch.float64)
        self.assertTrue(torch.autograd.gradcheck(
            lambda k,q:atom_scores(moment_state(k,values),q,codes,ridge=.1),
            (keys,queries),eps=1e-6,atol=1e-5,rtol=1e-4))

    def test_positive_appends_only_expand_compatibility_and_large_isotropic_mass_saturates(self):
        torch.manual_seed(95);codes=value_codes(dtype=torch.float64)
        keys=torch.rand(1,10,3,dtype=torch.float64);values=codes[torch.arange(10)%8][None]
        query=torch.rand(1,2,3,dtype=torch.float64)
        first=atom_scores(moment_state(keys[:,:4],values[:,:4]),query,codes)
        second=atom_scores(moment_state(keys,values),query,codes)
        self.assertTrue(torch.all(second<=first+1e-10))
        isotropic=(2*torch.eye(4,dtype=torch.float64))[None,None]
        self.assertTrue(torch.all(atom_scores(isotropic,torch.ones(1,1,1,dtype=torch.float64),codes)<=1))

    def test_identical_keys_with_conflicting_values_cannot_be_certified(self):
        codes=value_codes(dtype=torch.float64);keys=torch.ones(1,2,1,dtype=torch.float64)
        scores=atom_scores(moment_state(keys,codes[torch.tensor([[0,7]])]),keys[:,:1],codes)
        self.assertGreaterEqual(int((scores<=1+1e-10).sum()),2)

    def test_operator_is_permutation_invariant_for_additive_memory_and_differentiable(self):
        torch.manual_seed(93);operator=RetrievalOperator('atom').double()
        keys=torch.randn(2,5,16,dtype=torch.float64,requires_grad=True);queries=keys[:,:2];labels=torch.randint(8,(2,5))
        first,_=operator(keys,queries,labels);order=torch.tensor([4,1,0,3,2]);second,_=operator(keys[:,order],queries,labels[:,order])
        torch.testing.assert_close(first,second)
        first.sum().backward();self.assertTrue(torch.isfinite(keys.grad).all())
        self.assertTrue(torch.isfinite(operator.encoder.weight.grad).all())


if __name__=='__main__':unittest.main()
