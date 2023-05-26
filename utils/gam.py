from joblib import Parallel, delayed
from tqdm import tqdm

import numpy as np
from pygam import GAM


def _fit_gam(x, y):
    gam = GAM()
    gam.fit(X=x, y=y)
    
    return gam

def get_gams(adata, layer: str, time: str, n_jobs: int):
    cell_cycle_score = adata.obs[time].values
    cell_cycle_score = np.vstack([cell_cycle_score - 2 * np.pi, cell_cycle_score, cell_cycle_score + 2 * np.pi]).T
    
    res = Parallel(n_jobs)(
        delayed(_fit_gam)(
            cell_cycle_score, adata.layers[layer][:, var_id]
        ) for var_id in tqdm(range(adata.n_vars))
    )
    
    return dict(zip(adata.var_names, res))