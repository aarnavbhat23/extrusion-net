"""Separate source-to-result implementation, not external biological replication."""
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path
import re
from collections import defaultdict
import numpy as np
import pandas as pd

EX=Path(__file__).resolve().parent
report={}
trajectory=json.loads((EX/'results/trajectory_reversibility.json').read_text());checks=[]
for result in trajectory['rows']:
    source=EX/'external/trajectory_pilot'/result['file'];assert hashlib.sha256(source.read_bytes()).hexdigest()==result['sha256']
    tracks=defaultdict(list)
    with source.open() as handle:
        for r in csv.DictReader(handle,delimiter='\t'):
            date=re.match(r'^(\d{4}_\d{2}_\d{2})(?:_|$)',r['movie_index'])
            if date:tracks[(date[1],r['movie_index'],r['id'])].append((int(r['t']),float(r['dist'])))
    tracks={k:v for k,v in tracks.items() if len(set(t for t,_ in v))==len(v)}
    dates=sorted({k[0] for k in tracks},key=lambda d:hashlib.sha256(('reversibility-v1:'+d).encode()).hexdigest())
    dev=set(dates[:math.ceil(len(dates)/3)]);assert sorted(dev)==result['development_dates']
    thresholds=np.quantile([d for k,v in tracks.items() if k[0] in dev for _,d in v if np.isfinite(d)],[1/3,2/3])
    np.testing.assert_allclose(thresholds,result['development_tertiles_um'],rtol=0,atol=0)
    movie_values=defaultdict(list);saved={(r['date'],r['movie'],r['track']):r for r in result['tracks']}
    for key,values in tracks.items():
        if key[0] in dev:continue
        states={t:int(d>=thresholds[0])+int(d>=thresholds[1]) for t,d in values if np.isfinite(d)};scores=[]
        for index,lag in enumerate((1,2,4)):
            pairs=[(b,states[t+lag]) for t,b in states.items() if t+lag in states]
            total=sum((1 if (a,b) in ((0,1),(1,2),(2,0)) else -1 if a!=b else 0) for a,b in pairs)
            scores.append(total/len(pairs) if pairs else 0.)
            assert len(pairs)==saved[key]['pair_counts'][index]
        np.testing.assert_allclose(scores,saved[key]['scores'],atol=1e-15,rtol=0);movie_values[key[:2]].append(scores)
    dates_values=defaultdict(list)
    for (date,movie),values in movie_values.items():dates_values[date].append(np.mean(values,axis=0))
    x=np.array([np.mean(dates_values[date],axis=0) for date in sorted(dates_values)])
    np.testing.assert_allclose(x,[r['scores'] for r in result['groups']],atol=1e-15,rtol=0)
    norms=np.linalg.norm(x,axis=0);norms[norms==0]=1;stat=np.max(np.abs(x.sum(axis=0))/norms)
    # Complete reference enumeration; no production randomization function.
    signs=np.array(list(itertools.product((-1,1),repeat=len(x))))
    stats=np.max(np.abs(signs@x)/norms,axis=1);p=float(np.mean(stats>=stat-1e-12))
    assert p==result['test']['adjusted_p']
    checks.append({'file':result['file'],'confirmation_dates':len(x),'adjusted_p':p,'source_and_scores_verified':True})
report['extrusion']=checks

Path('results/trajectory_reversibility_verification.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
