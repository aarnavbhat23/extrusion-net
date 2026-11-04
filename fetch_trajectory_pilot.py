"""Verify two public trajectory file formats without fitting or scoring outcomes."""
import csv
import hashlib
import io
import json
from pathlib import Path
import ssl
import urllib.request
from collections import Counter

SOURCES=[{'record':5770531,'file':'Rad21_0_hr.tagged_set.tsv','delimiter':'\t','track':'id','frame':'t',
          'sha256':'7728e7ae178085cee27aea45c2247d673dff2f2749d256b537291e5dd4d0ca9b'},
         {'record':16949930,'file':'Adjacent_QC.csv','delimiter':';','track':'Track_pair','frame':'Frame',
          'sha256':'4ae83914fd272904b1933e30867be4d7e0a35a87353596d113565b61de56f923'}]

if __name__=='__main__':
    folder=Path('external/trajectory_pilot');folder.mkdir(parents=True,exist_ok=True);report=[]
    for source in SOURCES:
        url=f"https://zenodo.org/api/records/{source['record']}/files/{source['file']}/content"
        path=folder/source['file']
        if not path.exists():
            ca=Path('/etc/ssl/cert.pem');ctx=ssl.create_default_context(cafile=str(ca) if ca.exists() else None)
            path.write_bytes(urllib.request.urlopen(url,context=ctx,timeout=90).read())
        b=path.read_bytes()
        if hashlib.sha256(b).hexdigest()!=source['sha256']:raise ValueError('Source file checksum mismatch')
        reader=csv.DictReader(io.StringIO(b.decode()),delimiter=source['delimiter']);fields=reader.fieldnames;rows=list(reader)
        counts=Counter(r[source['track']] for r in rows)
        report.append({**source,'url':url,'license':'CC BY 4.0','bytes':len(b),'columns':fields,
            'rows':len(rows),'tracks':len(counts),'minimum_rows_per_track':min(counts.values()),
            'maximum_rows_per_track':max(counts.values()),
            'movies':len({r['movie_index'] for r in rows}) if 'movie_index' in fields else None,
            'scope':'Metadata-only feasibility; no outcome distributions fitted or evaluated.'})
    Path('results/trajectory_feasibility.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
