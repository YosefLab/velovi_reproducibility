"""
Code taken and adapted from `https://github.com/pachterlab/GFCP_2022`

BSD 2-Clause License

Copyright (c) 2020, Pachter Lab
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


import numpy as np
import scipy
from scipy.fft import irfft2
from scipy.stats import rv_discrete


def simulate_burst_model(nCells=2000, nGenes=100, T=20, n_cell_types=10, seed=42):
    X = np.zeros((2, nCells, nGenes))

    n_cells_per_type = nCells // n_cell_types
    g_true = np.zeros((nGenes, n_cell_types))  # deg rate
    b_true = np.zeros(nGenes)  # burst size

    cell_types = np.zeros(nCells)
    np.random.seed(seed)
    K = np.zeros((nGenes, n_cell_types))
    for j in range(nGenes):
        beta = 1
        b = np.clip(10 ** np.random.normal(0.3, 0.8), 0.05, 25)
        b_true[j] = b
        gamma_mean = np.random.normal(-0.3, 0.3)

        for i in range(n_cell_types):
            FLAG = True

            while FLAG:
                kini = np.clip(10 ** np.random.normal(-1, 0.5), 0.005, 1)
                K[j, i] = kini
                gamma = np.clip(10 ** np.random.normal(gamma_mean, 0.1), 0.08, 4)
                g_true[j, i] = gamma

                if (i + 1) < n_cell_types:
                    IND = range(i * n_cells_per_type, (i + 1) * n_cells_per_type)
                else:
                    IND = range(i * n_cells_per_type, nCells)
                mu = np.array([kini * b / beta, kini * b / gamma])
                std = np.sqrt(mu + mu * [b, b * beta / (beta + gamma)])
                lm = np.clip(np.ceil(mu + 5 * std), 10, np.inf)
                lm = [int(x) for x in lm]

                if np.prod(lm) > 15000:
                    pass
                else:
                    mesh = np.meshgrid(*[np.arange(x) for x in lm], indexing="ij")
                    mesh = (mesh[0].flatten(), mesh[1].flatten())
                    var_val = np.asarray(
                        [(mesh[0][i], mesh[1][i]) for i in range(np.prod(lm))]
                    )
                    prob = cme_integrator(np.array([kini, b, beta, gamma]), lm)
                    prob = prob.flatten()
                    var_inds = np.arange(0, len(var_val))
                    distrib = rv_discrete(values=(var_inds, prob))
                    VAL = var_val[distrib.rvs(size=len(IND))]
                    X[:, IND, j] = VAL.T
                    FLAG = False

            if j == 0:
                cell_types[IND] = i
    return X, cell_types, cell_types, K, g_true[:, 0], b_true


def cme_integrator(
    p, lm, method="fixed_quad", fixed_quad_T=10, quad_order=60, quad_vec_T=np.inf
):
    kini, b, bet, gam = p
    u = []
    mx = np.copy(lm)

    # initialize the generating function evaluation points
    mx[-1] = mx[-1] // 2 + 1
    for i in range(len(mx)):
        l = np.arange(mx[i])
        u_ = np.exp(-2j * np.pi * l / lm[i]) - 1
        u.append(u_)
    g = np.meshgrid(*[u_ for u_ in u], indexing="ij")
    for i in range(len(mx)):
        g[i] = g[i].flatten()[:, np.newaxis]

    if bet != gam:  # compute weights for the ODE solution.
        f = b * bet / (bet - gam)
        g[1] *= f
        g[0] *= b
        g[0] -= g[1]
    else:
        g[1] *= b * gam
        g[0] *= b

    # define function to integrate by quadrature.
    fun = lambda x: INTFUN(x, g, bet, gam)
    if method == "quad_vec":
        T = quad_vec_T * (1 / bet + 1 / gam + 1 / kini)
        gf = scipy.integrate.quad_vec(fun, 0, T)[0]
    if method == "fixed_quad":
        T = fixed_quad_T * (1 / bet + 1 / gam + 1 / kini)
        gf = scipy.integrate.fixed_quad(fun, 0, T, n=quad_order)[0]

    # convert back to the probability domain, renormalize to ensure non-negativity.
    gf = np.exp(
        kini * gf
    )  # gf can be multiplied by k in the argument, but this is not relevant for the 3-parameter input.
    gf = gf.reshape(tuple(mx))
    Pss = irfft2(gf, s=tuple(lm))
    Pss[Pss < 1e-9] = 1e-9
    Pss = np.abs(Pss) / np.sum(np.abs(Pss))
    return Pss


def INTFUN(x, g, bet, gam):
    U = np.exp(-bet * x) * g[0] + np.exp(-gam * x) * g[1]
    return U / (1 - U)
