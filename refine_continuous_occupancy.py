"""Prespecified identity-based reanalysis; never overwrite original interval results."""
import json,time
from pathlib import Path
import numpy as np
from continuous_occupancy import enclose

ROOT=Path(__file__).resolve().parent


def main():
    source=json.loads((ROOT/'results/continuous_occupancy_study.json').read_text());rows=[];arrays={}
    for index,row in enumerate(source['rows']):
        observations=[(x['gamma'],np.array(x['counts'])) for x in row['observations']]
        start=time.perf_counter();result,lo,hi=enclose(observations,aggregate=True);seconds=time.perf_counter()-start
        truth=np.array(row['exact_target']);covered=not result['rejected'] and bool(np.all((truth>=np.array(result['lower'])-1e-12)&(truth<=np.array(result['upper'])+1e-12)))
        contains=bool(np.any(np.all((np.array(row['parameters'])>=lo)&(np.array(row['parameters'])<=hi),axis=1)))
        rows.append({'parameters':row['parameters'],'repeat':row['repeat'],'continuous':result,'seconds':seconds,
                     'covers_generating_target':covered,'contains_generating_parameters':contains,
                     'original_max_width':row['continuous']['max_target_width']})
        arrays[f'case_{index}_low']=lo;arrays[f'case_{index}_high']=hi
        (ROOT/'results/continuous_occupancy_refinement.json').write_text(json.dumps({'amendment_commit':'d1cf3e6','scope':'Development reanalysis of the same six cases; not fresh confirmation.','rows':rows},indent=2)+'\n')
        np.savez_compressed(ROOT/'results/continuous_occupancy_refined_boxes.npz',**arrays)
        print(index,'covered',covered,'parameters',contains,'attained',result['attained'],'width',result.get('max_target_width'),
              'boxes',result.get('retained_boxes'),'seconds',seconds,flush=True)


if __name__=='__main__':main()
