import numpy as np
import pandas as pd
import time
from PIL import Image
from scipy.sparse.linalg import svds
from scipy.spatial.distance import pdist, squareform


def pixel_intensity(adata, key=None, windowsize=20):
    t1 = time.time()
    if key is None:
        key = list(adata.uns['spatial'].keys())[0]
    scale = adata.uns['spatial'][key]['scalefactors']['tissue_hires_scalef']

    image = np.uint8(adata.uns['spatial'][key]["images"]['hires']*255)
    if image.shape[-1]>3:
        image = image[:,:,:3]

    image = Image.fromarray(image)
    image.convert("L")
    intensity_3d = np.zeros((len(adata.obs_names), 3))
    intensity_1d = np.zeros(len(adata.obs_names))
    x = np.round(adata.obsm['spatial'][:,1]*scale)
    y = np.round(adata.obsm['spatial'][:,0]*scale)
    for i in (range(len(adata.obs_names))):
        subimage = np.asarray(image.crop((y[i]-windowsize, x[i]-windowsize, y[i]+windowsize, x[i]+windowsize)))
        intensity_1d[i] = subimage.mean()
        intensity_3d[i, :] = subimage.mean(axis=0).mean(axis=0).T
    adata.obs['pixel'] = intensity_1d
    adata.obsm['pixel'] = ( intensity_3d - intensity_3d.mean(axis=0) ) / intensity_3d.std(axis=0)
    t2 = time.time()
    print('Time elapsed: %.2f seconds'%(t2-t1))

def make_kernel(adata, n=100, bandwidth=1, im_feats_weight=0.3):
    t1 = time.time()
    X = np.concatenate((adata.obsm['spatial'][:, 0:2], adata.obsm['pixel']), axis=1)
    x = X.T
    x = (x-x.mean(axis=1).reshape(-1,1))/x.std(axis=1).reshape(-1,1)
    Xn = x.T
    if X.shape[1] > 2:
        Xn[:, 2:X.shape[1]] = im_feats_weight*Xn[:, 2:X.shape[1]]
    pw_dist = squareform(pdist(Xn, 'euclidean'))
    adata.obsp['pw_dist'] = pw_dist
    # adata.obsp['kernel'] = (1/(np.sqrt(2)*np.pi*bandwidth)**X.shape[1]) * np.exp(-pw_dist**2 / (2 * bandwidth**2))
    adata.obsp['kernel'] = (1/(np.sqrt(2*np.pi)*bandwidth)**X.shape[1]) * np.exp(-pw_dist**2 / (2 * bandwidth**2))
    while n>Xn.shape[0]:
        n -= 50
    u,s,v = svds(adata.obsp['kernel'], n)
    adata.obsm['kernel'] = u.dot(np.diag(s))
    t2 = time.time()
    print('Time elapsed: %.2f seconds'%(t2-t1))