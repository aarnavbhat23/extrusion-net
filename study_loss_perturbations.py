"""Structural calibration screen with regulated removal and unaffected loss."""
import argparse
import json
from pathlib import Path
import numpy as np
from clock_calibration import flow_snapshot, generator_snapshot

ROOT = Path(__file__).resolve().parent
MECHANISMS = ((1.15,.85,.27,.35,.65,.45,.15),
              (.75,1.25,.5,.6,.3,.9,.15), (1.3,.7,.7,.45,.55,.6,.15))
# Conditions are (relative speed, turnover multiplier, observation rate).
DESIGNS = {
    'relative_speed': ([(s,1.,.3) for s in (.5,1.,2.)],True,True),
    'turnover_known_zero_loss': ([(1.,f,.3) for f in (.5,1.,2.)],True,False),
    'turnover_unknown_loss': ([(1.,f,.3) for f in (.5,1.,2.)],True,True),
    'two_clocks_unknown_loss': ([(1.,1.,g) for g in (.1,.3)],True,True),
    'three_conditions': ([(1.,1.,.1),(1.,1.,.3),(1.,2.,.3)],True,True),
    'absolute_two_turnovers': ([(1.,1.,.3),(1.,2.,.3)],False,True),
}


def snapshot(parameters, condition, conditional=True, independent=False):
    vl,vr,mu,pl,pr,h,delta = parameters
    speed,factor,gamma = condition
    removal = factor*mu+delta
    effective = (speed*vl,speed*vr,removal,pl,pr,(factor*mu*h+delta)/removal)
    p = generator_snapshot(effective,gamma) if independent else flow_snapshot(effective,gamma)
    return p[:-1]/p[:-1].sum() if conditional else p


def effective_rates(q, gamma, vl, vr):
    p = np.asarray(q)[:81].reshape(9,9)
    a,b = p[1,0]/p[0,0],p[0,1]/p[0,0]
    removal = vl/a-vl-vr-gamma
    pl = 1-p[4,0]/(p[2,0]*a*a)
    pr = 1-p[0,4]/(p[0,2]*b*b)
    captured = None
    if pl > 1e-10:
        d = p[3,0]/p[2,0]-(1-pl)*a
        captured = vl*pl/d-vr-gamma
    elif pr > 1e-10:
        d = p[0,3]/p[0,2]-(1-pr)*b
        captured = vr*pr/d-vl-gamma
    return removal,pl,pr,captured


def recover(observations, absolute=False):
    """Exact-probability inverse for the two specified sufficient designs only."""
    if absolute:
        q1,q2 = observations
        rate = .3/q1[0]
        vl,vr = rate*q1[9]/q1[0],rate*q1[1]/q1[0]
        first,second = effective_rates(q1,.3,vl,vr),effective_rates(q2,.3,vl,vr)
    else:
        q1,q2,q3 = observations
        vl = .2/(q2[0]/q2[9]-q1[0]/q1[9])
        vr = vl*q1[1]/q1[9]
        first,second = effective_rates(q2,.3,vl,vr),effective_rates(q3,.3,vl,vr)
    mu = second[0]-first[0]
    delta = first[0]-mu
    h = (first[3]-delta)/mu if first[3] is not None else None
    return [float(vl),float(vr),float(mu),float(first[1]),float(first[2]),
            float(h) if h is not None else None,float(delta)]


def sensitivity(parameters, conditions, conditional, unknown_loss, step):
    columns = []
    for j in range(7 if unknown_loss else 6):
        above,below = np.array(parameters),np.array(parameters)
        above[j] *= 1+step; below[j] *= 1-step
        plus = np.concatenate([snapshot(above,c,conditional) for c in conditions])
        minus = np.concatenate([snapshot(below,c,conditional) for c in conditions])
        columns.append((plus-minus)/(2*step))
    singular = np.linalg.svd(np.array(columns).T,compute_uv=False)
    return {'relative_step':step,'singular_values':singular.tolist(),
            'diagnostic_rank':int(np.sum(singular > singular[0]*1e-7)),
            'parameter_count':len(columns)}


def run():
    rows = []; witnesses = []; inverse_errors = []
    for case, base in enumerate(MECHANISMS):
        for name,(conditions,conditional,unknown_loss) in DESIGNS.items():
            p = base if unknown_loss else (*base[:6],0.)
            observations = [snapshot(p,c,conditional,True) for c in conditions]
            for c,q in zip(conditions,observations):
                np.testing.assert_allclose(q,snapshot(p,c,conditional),atol=1e-12,rtol=0)
            row = {'case':case,'design':name,'parameters':p,'conditions':conditions,
                   'conditional':conditional,
                   'sensitivities':[sensitivity(p,conditions,conditional,unknown_loss,h) for h in (1e-4,5e-5)]}
            if name in ('three_conditions','absolute_two_turnovers'):
                reconstructed = recover(observations,name == 'absolute_two_turnovers')
                error = float(np.max(np.abs(np.array(reconstructed)-p)))
                row.update(recovered=reconstructed,maximum_recovery_error=error)
                inverse_errors.append(error)
            rows.append(row)
        vl,vr,mu,pl,pr,h,delta = base; scale=1.6
        # Uniform background loss absorbs the offset for every speed/turnover input.
        transformed = (scale*vl,scale*vr,scale*mu,pl,pr,h,scale*delta+(scale-1)*.3)
        for name in ('relative_speed','turnover_unknown_loss'):
            distances = [float(np.abs(snapshot(base,c,True,True)-snapshot(transformed,c,True,True)).sum()/2)
                         for c in DESIGNS[name][0]]
            assert max(distances) < 1e-12
            witnesses.append({'case':case,'kind':name,'first':base,'second':transformed,'condition_total_variations':distances})
        # Same effective losses at f=1, at every clock; different regulated response.
        shift=.05
        decomposition = (vl,vr,mu-shift,pl,pr,(mu*h-shift)/(mu-shift),delta+shift)
        distances = [float(np.abs(snapshot(base,c,True,True)-snapshot(decomposition,c,True,True)).sum()/2)
                     for c in DESIGNS['two_clocks_unknown_loss'][0]]
        target_distance = float(np.abs(snapshot(base,(1.,2.,.3),True,True)-snapshot(decomposition,(1.,2.,.3),True,True)).sum()/2)
        assert max(distances) < 1e-12 and target_distance > 0
        witnesses.append({'case':case,'kind':'loss_decomposition','first':base,'second':decomposition,
                          'condition_total_variations':distances,'turnover_target_total_variation':target_distance})
    return {'protocol_commit':'062bef1','scope':'Structural matched-family feasibility, no finite-data or biological validation.',
            'rows':rows,'explicit_witnesses':witnesses,
            'maximum_inverse_error':max(inverse_errors),
            'structural_continuation_gate_passed':max(inverse_errors) < 1e-9}


if __name__ == '__main__':
    parser=argparse.ArgumentParser(); parser.add_argument('--verify',action='store_true'); args=parser.parse_args()
    result=run(); path=ROOT/'results/loss_perturbations.json'
    if args.verify:
        saved=json.loads(path.read_text())
        assert saved['protocol_commit'] == result['protocol_commit']
        assert len(saved['rows']) == len(result['rows'])
        for a,b in zip(saved['rows'],result['rows']):
            assert (a['case'],a['design']) == (b['case'],b['design'])
            np.testing.assert_allclose(a['parameters'],b['parameters'],atol=0,rtol=0)
            np.testing.assert_allclose(a['conditions'],b['conditions'],atol=0,rtol=0)
            for x,y in zip(a['sensitivities'],b['sensitivities']):
                assert x['diagnostic_rank'] == y['diagnostic_rank']
                np.testing.assert_allclose(x['singular_values'],y['singular_values'],atol=1e-8,rtol=1e-5)
            if 'recovered' in a: np.testing.assert_allclose(a['recovered'],b['recovered'],atol=1e-9,rtol=0)
        assert len(saved['explicit_witnesses']) == len(result['explicit_witnesses']) == 9
        for a,b in zip(saved['explicit_witnesses'],result['explicit_witnesses']):
            assert (a['case'],a['kind']) == (b['case'],b['kind'])
            np.testing.assert_allclose(a['first'],b['first'],atol=0,rtol=0)
            np.testing.assert_allclose(a['second'],b['second'],atol=1e-14,rtol=0)
            np.testing.assert_allclose(a['condition_total_variations'],b['condition_total_variations'],atol=1e-12,rtol=0)
            if 'turnover_target_total_variation' in a:
                np.testing.assert_allclose(a['turnover_target_total_variation'],b['turnover_target_total_variation'],atol=1e-12,rtol=0)
        np.testing.assert_allclose(saved['maximum_inverse_error'],result['maximum_inverse_error'],atol=1e-12,rtol=0)
        assert saved['structural_continuation_gate_passed'] == result['structural_continuation_gate_passed']
    else:
        path.write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print('Maximum inverse error:',result['maximum_inverse_error'])
    for name in DESIGNS:
        print(name,[(r['case'],[s['diagnostic_rank'] for s in r['sensitivities']]) for r in result['rows'] if r['design']==name])
