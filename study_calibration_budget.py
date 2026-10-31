"""Local Gaussian information bounds for balanced extrusion routing."""
import json
from pathlib import Path

import numpy as np

ROOT=Path(__file__).resolve().parent


def tangent_contacts(gates, pass_rate):
    """Jacobian along l=a+z, r=a-z, q=1 at z=0; all upper contacts."""
    size=gates+1
    rows=[]
    for gap in range(2,size):
        coefficients=pass_rate**(gap-1)*(gap-1-2*np.arange(gap))
        for start in range(size-gap):
            row=np.zeros(gates)
            row[start:start+gap]=coefficients
            rows.append(row)
    return np.asarray(rows).reshape(-1,gates)


def laplacian(gates):
    differences=np.diff(np.eye(gates),axis=0)
    return differences.T@differences


def diffusion_bound(pass_rate):
    r=pass_rate**2
    return r*(1+r)*(1+8*r+r*r)/(1-r)**7


def information_covariance(jacobian,sites):
    info=jacobian.T@jacobian
    info[np.asarray(sites),np.asarray(sites)]+=1
    return np.linalg.solve(info,np.eye(len(info)))


def main():
    rows=[]
    for pass_rate in (.3,.6,.8):
        for gates in (8,16,32,64,128,256):
            jac=tangent_contacts(gates,pass_rate);fisher=jac.T@jac
            path=laplacian(gates);constant=diffusion_bound(pass_rate)
            lower_eig=float(np.linalg.eigvalsh(fisher-pass_rate**2*path).min())
            upper_eig=float(np.linalg.eigvalsh(constant*path-fisher).min())
            for budget in sorted(set([1,2,4,8,max(1,gates//8)])):
                if budget>gates:continue
                for placement in ('boundary_cluster','evenly_spaced'):
                    sites=(np.arange(budget) if placement=='boundary_cluster' else
                           np.floor((np.arange(budget)+.5)*gates/budget).astype(int))
                    covariance=information_covariance(jac,sites)
                    row={'gates':gates,'pass_rate':pass_rate,'budget':budget,'placement':placement,
                         'sites':sites.tolist(),'mean_variance_over_sigma_squared':float(np.trace(covariance)/gates),
                         'max_variance_over_sigma_squared':float(np.diag(covariance).max()),
                         'any_placement_lower_bound':float((gates-budget)**2/(6*constant*gates*(budget+1))),
                         'spectral_lower_min_eigenvalue':lower_eig,'spectral_upper_min_eigenvalue':upper_eig}
                    if budget==1 and placement=='boundary_cluster':
                        row['boundary_lower_bound']=1+(gates-1)/(2*constant)
                        row['boundary_upper_bound']=1+(gates-1)/(2*pass_rate**2)
                    rows.append(row)
            print('pass',pass_rate,'gates',gates,flush=True)
    # Independent observation-space simulation in the linearized Gaussian model.
    # This checks the computation, not nonlinear or neural-network performance.
    rng=np.random.default_rng(20261101);checks=[]
    for budget in (1,4,8):
        gates=32;jac=tangent_contacts(gates,.6)
        sites=np.floor((np.arange(budget)+.5)*gates/budget).astype(int)
        design=np.vstack([jac,np.eye(gates)[sites]])
        covariance=information_covariance(jac,sites)
        errors=np.linalg.solve(design.T@design,design.T@rng.normal(size=(len(design),2000)))
        checks.append({'gates':gates,'budget':budget,'repetitions':2000,
                       'predicted_mean_mse':float(np.trace(covariance)/gates),
                       'observed_mean_mse':float(np.mean(errors**2))})
    result={'scope':'Local Cramer-Rao analysis at homogeneous interior gates, q=1 and known l+r; independent homoscedastic Gaussian contact and calibration noise. Not a global minimax or learned-network guarantee.',
            'rows':rows,'linearized_simulation':checks}
    (ROOT/'results/calibration_budget_study.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
