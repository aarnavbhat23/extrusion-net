"""Frozen real-trajectory screen, with outcomes clustered by acquisition date."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import ssl
import urllib.request
from collections import defaultdict
import numpy as np
from trajectory_reversibility import track_scores,reversal_test

SOURCES=[('Rad21_0_hr.unfiltered.tagged_set.tsv','7e9a1063b1035b5c380a0a45c067a311'),
         ('Rad21_0_hr.tagged_set.tsv','fd306f4feb651b823ea41a328c331746')]


def load(path):
    tracks=defaultdict(list);excluded=set();rows=0
    with path.open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            rows+=1;movie=row['movie_index'];match=re.match(r'^(\d{4}_\d{2}_\d{2})(?:_|$)',movie)
            if not match:excluded.add(movie);continue
            tracks[(match[1],movie,row['id'])].append((float(row['t']),float(row['dist'])))
    ambiguous=[{'date':key[0],'movie':key[1],'track':key[2],'rows':len(values)} for key,values in sorted(tracks.items())
               if len({t for t,_ in values})!=len(values)]
    for item in ambiguous:del tracks[(item['date'],item['movie'],item['track'])]
    return tracks,{'rows':rows,'excluded_unmapped_movies':sorted(excluded),'excluded_duplicate_tracks':ambiguous}


def analyze(path):
    tracks,metadata=load(path);dates=sorted({key[0] for key in tracks},
        key=lambda date:hashlib.sha256(('reversibility-v1:'+date).encode()).hexdigest())
    dev=set(dates[:math.ceil(len(dates)/3)]);test=set(dates)-dev
    distances=[distance for key,values in tracks.items() if key[0] in dev for _,distance in values if np.isfinite(distance)]
    thresholds=np.quantile(distances,[1/3,2/3]);movies=defaultdict(list);track_rows=[]
    for (date,movie,track),values in sorted(tracks.items()):
        frames,dist=np.array(values).T;scores,counts=track_scores(frames,dist,thresholds)
        # Independent path reversal check includes missing frames and unequal track lengths.
        reverse,_=track_scores(frames.min()+frames.max()-frames,dist,thresholds)
        if not np.allclose(reverse,-scores,atol=1e-14):raise AssertionError('Score is not antisymmetric')
        if date in test:
            movies[(date,movie)].append(scores)
            track_rows.append({'date':date,'movie':movie,'track':track,'scores':scores.tolist(),'pair_counts':counts})
    groups=defaultdict(list)
    for (date,movie),values in sorted(movies.items()):groups[date].append(np.mean(values,axis=0))
    group_rows=[{'date':date,'scores':np.mean(values,axis=0).tolist(),'movies':len(values)} for date,values in sorted(groups.items())]
    result=reversal_test([r['scores'] for r in group_rows])
    return {**metadata,'development_dates':sorted(dev),'confirmation_dates':sorted(test),
            'development_tertiles_um':thresholds.tolist(),'confirmation_movies':len(movies),
            'confirmation_tracks':len(track_rows),'groups':group_rows,'tracks':track_rows,'test':result}


if __name__=='__main__':
    folder=Path('external/trajectory_pilot');folder.mkdir(parents=True,exist_ok=True)
    report={'protocol_commit':'ac2578f','input_amendment':'docs/TRAJECTORY_INPUT_AMENDMENT.md','cadence_seconds':20,'lags_frames':[1,2,4],
            'source_doi':'10.5281/zenodo.5770531','license':'CC BY 4.0','rows':[]}
    for name,md5 in SOURCES:
        path=folder/name;url=f'https://zenodo.org/api/records/5770531/files/{name}/content'
        if not path.exists():
            ca=Path('/etc/ssl/cert.pem');ctx=ssl.create_default_context(cafile=str(ca) if ca.exists() else None)
            path.write_bytes(urllib.request.urlopen(url,context=ctx,timeout=90).read())
        blob=path.read_bytes()
        if hashlib.md5(blob).hexdigest()!=md5:raise ValueError('Source MD5 mismatch')
        row={'file':name,'url':url,'sha256':hashlib.sha256(blob).hexdigest(),'md5':md5,
             'primary':'unfiltered' in name,**analyze(path)};report['rows'].append(row)
        print(name,json.dumps(row['test']),flush=True)
    Path('results/trajectory_reversibility.json').write_text(json.dumps(report,indent=2)+'\n')
