"""Post-screen explanation, frozen trained encoders; not new confirmation."""
import json
from pathlib import Path
import numpy as np
import torch
from moment_memory import RetrievalOperator,moment_state,atom_scores
from study_moment_retrieval import pool,episode

ROOT=Path(__file__).resolve().parent


@torch.no_grad()
def main():
    torch.set_num_threads(2)
    keys,labels=episode(pool(9202,1024),np.random.default_rng(9660),64,64)
    report={'scope':'Exploratory explanation after failed retrieval gate; query first eight keys of nested prefixes.',
            'seeds':[0,1,2],'rows':[]}
    for seed in report['seeds']:
        model=RetrievalOperator('atom')
        model.load_state_dict(torch.load(ROOT/'results'/f'moment_operator_seed{seed}.pt',weights_only=True)['state_dict'])
        model.double();features=model.features(keys.double());queries=features[:,:8];prior=None
        for n in (8,16,32,64):
            state=moment_state(features[:,:n],model.codes[labels[:,:n]])
            scores=atom_scores(state,queries,model.codes)
            gram=torch.einsum('bqd,bdij->bqij',queries,state)+1e-3*torch.eye(4,dtype=torch.float64)
            minimum=torch.linalg.eigvalsh(gram)[...,0]
            row={'seed':seed,'bindings':n,'mean_compatible_values':float((scores<=1+1e-8).sum(-1).double().mean()),
                 'all_values_compatible_fraction':float((scores<=1+1e-8).all(-1).double().mean()),
                 'lambda_min_at_least_two_fraction':float((minimum>=2).double().mean()),
                 'largest_score_increase_from_previous':None if prior is None else float((scores-prior).max())}
            if prior is not None:assert torch.all(scores<=prior+1e-8)
            report['rows'].append(row);prior=scores
    (ROOT/'results/moment_saturation_study.json').write_text(json.dumps(report,indent=2)+'\n')
    for n in (8,16,32,64):
        rows=[r for r in report['rows'] if r['bindings']==n]
        print(n,{k:np.mean([r[k] for r in rows]) for k in ('mean_compatible_values','all_values_compatible_fraction','lambda_min_at_least_two_fraction')})


if __name__=='__main__':main()
