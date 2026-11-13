"""Independent decimal flow verifies saved feasible pairs, without optimizing."""
from collections import defaultdict
from decimal import Decimal,localcontext
import json
from pathlib import Path
import numpy as np
from continuous_occupancy import calibration_intervals

ROOT=Path(__file__).resolve().parent


def decimal_snapshot(parameters,gamma=.3,factor=1.):
    d=lambda x:Decimal.from_float(float(x))
    mu,pl,pr,h=map(d,parameters);mu*=d(factor);gamma=d(gamma)
    arrivals=defaultdict(Decimal);arrivals[0,0,0,0]=Decimal(1);result=[Decimal(0)]*81
    for a in range(9):
        for b in range(9):
            for cl in range(2):
                for cr in range(2):
                    incoming=arrivals[a,b,cl,cr]
                    if not incoming:continue
                    left=int(a<8 and not cl);right=int(b<8 and not cr)
                    residence=incoming/(left+right+gamma+mu*(h if cl or cr else 1))
                    result[a*9+b]+=gamma*residence
                    if left:
                        p=pl if a+1==3 else Decimal(0)
                        arrivals[a+1,b,0,cr]+=residence*(1-p);arrivals[a+1,b,1,cr]+=residence*p
                    if right:
                        p=pr if b+1==3 else Decimal(0)
                        arrivals[a,b+1,cl,0]+=residence*(1-p);arrivals[a,b+1,cl,1]+=residence*p
    return result+[Decimal(1)-sum(result)]


def main():
    source=json.loads((ROOT/'results/continuous_occupancy_study.json').read_text())
    diagnostic=json.loads((ROOT/'results/occupancy_ambiguity_diagnostic.json').read_text());rows=[]
    with localcontext() as context:
        context.prec=70
        for data,row in zip(source['rows'],diagnostic['rows']):
            targets=[];slacks=[]
            for key in ('lower_absence_witness','upper_absence_witness'):
                parameters=row[key]['parameters'];minimum=Decimal(1)
                for observation in data['observations']:
                    lo,hi=calibration_intervals(np.array(observation['counts']))
                    predicted=decimal_snapshot(parameters,observation['gamma'])
                    minimum=min(minimum,min(p-Decimal.from_float(float(l)) for p,l in zip(predicted,lo)),
                                min(Decimal.from_float(float(u))-p for p,u in zip(predicted,hi)))
                if minimum<0:raise AssertionError('Witness violates represented calibration constraints')
                slacks.append(str(minimum));targets.append(decimal_snapshot(parameters,factor=2.)[-1])
            span=targets[1]-targets[0]
            if span<=Decimal('.02'):raise AssertionError('Witness span does not obstruct stated precision')
            rows.append({'parameters':data['parameters'],'repeat':data['repeat'],'minimum_calibration_slacks':slacks,
                         'attained_absence_span':str(span),'exceeds_point02':True})
    (ROOT/'results/occupancy_witness_verification.json').write_text(json.dumps({'decimal_digits':70,'scope':'Verification of supplied numerical witnesses against represented rectangular calibration bounds; no global optimum claim.','rows':rows},indent=2)+'\n')
    print('All 12 witness populations satisfy every calibration inequality at 70-digit precision; all six spans exceed .02.')


if __name__=='__main__':main()
