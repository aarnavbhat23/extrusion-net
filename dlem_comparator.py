"""Pinned authors' forward for a common random-delay observable; no redistribution."""
import hashlib,importlib.util,json,ssl,urllib.request
from pathlib import Path
import numpy as np
import torch

ROOT=Path(__file__).resolve().parent
REVISION='e55f9226192c34c34990518713a15ec42dff884a'
URL=f'https://raw.githubusercontent.com/chikinalab/dLEMpytorch/{REVISION}/dlem/models/netdlem2.py'
SHA256='35fcb75676f895d22e801520be9ae21001a172b66bcbbc2ef36ccfe17e94d290'


def fetch():
    path=ROOT/'external/netdlem2.py';path.parent.mkdir(exist_ok=True)
    if not path.exists():
        data=urllib.request.urlopen(URL,context=ssl.create_default_context(cafile='/etc/ssl/cert.pem')).read()
        if hashlib.sha256(data).hexdigest()!=SHA256:raise ValueError('Published source hash mismatch')
        path.write_bytes(data)
    return path


def load_model():
    path=ROOT/'external/netdlem2.py'
    if not path.exists():raise FileNotFoundError('Run python dlem_comparator.py before the comparator test')
    if hashlib.sha256(path.read_bytes()).hexdigest()!=SHA256:raise ValueError('Published source changed')
    spec=importlib.util.spec_from_file_location('published_netdlem2',path);module=importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.DLEM(17).double()


def snapshot(parameters,gamma=.3,removal_factor=1.,model=None):
    model=load_model() if model is None else model
    mu=parameters[16]*removal_factor;left=np.zeros(17);right=np.zeros(17)
    left[8-np.arange(8)]=parameters[:8];right[8+np.arange(8)]=parameters[8:16]
    scale=1/max(float(max(left.max(),right.max(),mu+gamma)),1.)
    # In this pinned PyTorch version left drives the right anchor, right the left.
    with torch.no_grad():
        model.left.copy_(torch.tensor(right*scale));model.right.copy_(torch.tensor(left*scale))
        model.unload.fill_(1-np.sqrt((mu+gamma)*scale));model.const.fill_(1.)
        initial=torch.zeros((1,17),dtype=torch.float64);initial[0,8]=gamma/(left[8]+right[8]+mu+gamma)
        matrix=model.contact_map_prediction(initial)[0].numpy()
    grid=matrix[np.ix_(8-np.arange(9),8+np.arange(9))]
    return np.r_[grid.ravel(),max(0.,1-grid.sum())]


if __name__=='__main__':
    path=fetch()
    print(json.dumps({'path':str(path),'source':URL,'revision':REVISION,'sha256':SHA256,
                      'attribution':'Chikina lab and dLEM authors; unchanged source, not redistributed.'},indent=2))
