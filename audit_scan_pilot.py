"""Replay every pilot checkpoint and retain token-level predictions for inspection."""
import hashlib
import json
from pathlib import Path

import numpy as np
import torch

import task
from models import SequenceModel
from train import TEST_LENGTHS


def main():
    torch.set_num_threads(2)
    directory=Path('results/scan_pilot')
    rows=[json.loads(line) for line in (directory/'runs.jsonl').read_text().splitlines()]
    audit=[]
    for row in rows:
        name=f"{row['mixer']}_s{row['seed']}"
        model=SequenceModel(row['mixer'],layers=1)
        model.load_state_dict(torch.load(directory/f'{name}.pt',weights_only=True))
        model.eval()
        predictions={}
        checks={}
        for n,(seed,count) in TEST_LENGTHS.items():
            tokens,labels=task.make_split(seed,count,[n])[n]
            with torch.no_grad():
                predicted=np.concatenate([model(torch.from_numpy(tokens[start:start+16])).argmax(-1).numpy()
                                          for start in range(0,len(tokens),16)])
            mask=labels!=task.IGNORE
            accuracy=float((predicted[mask]==labels[mask]).mean())
            assert round(accuracy,4)==row['test_acc'][str(n)],(name,n,accuracy)
            predictions[f'pred_{n}']=predicted.astype(np.int16)
            predictions[f'target_{n}']=labels.astype(np.int16)
            checks[str(n)]={'correct':int((predicted[mask]==labels[mask]).sum()),'total':int(mask.sum()),'accuracy':accuracy}
        destination=directory/f'{name}_predictions.npz'
        np.savez_compressed(destination,**predictions)
        audit.append({'run':name,'checks':checks,'prediction_sha256':hashlib.file_digest(destination.open('rb'),'sha256').hexdigest()})
    (directory/'replay_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
    print('All nine checkpoints reproduced at all four test lengths.')


if __name__=='__main__':
    main()
