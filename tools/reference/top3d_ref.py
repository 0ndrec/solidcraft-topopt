"""Independent reference for topopt-core: assembled sparse stiffness, direct solve (scipy),
the same hat density filter and OC update, no projection. Written separately from the Rust code
(its own Hex8 stiffness and node numbering) to cross-check the matrix-free + multigrid solver,
the filter and the optimality-criteria loop.

    tools/reference/.venv/Scripts/python tools/reference/top3d_ref.py NX NY NZ ITERS [volfrac] [rmin]

Problem (top3d cantilever): x = 0 face clamped, total load -1 (z) spread evenly over the nodes of
the edge x = L, z = 0; E = 1, nu = 0.3, element size 1. Prints one JSON line per iteration.
"""
import json
import sys

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import spsolve


def hex8_ke(nu):
    lam = nu / ((1 + nu) * (1 - 2 * nu))
    mu = 1 / (2 * (1 + nu))
    D = np.full((3, 3), lam) + 2 * mu * np.eye(3)
    D = np.block([[D, np.zeros((3, 3))], [np.zeros((3, 3)), mu * np.eye(3)]])
    corners = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0], [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]])
    s = 2 * corners - 1
    KE = np.zeros((24, 24))
    g = 1 / np.sqrt(3)
    for gx in (-g, g):
        for gy in (-g, g):
            for gz in (-g, g):
                dN = np.stack([
                    s[:, 0] * (1 + s[:, 1] * gy) * (1 + s[:, 2] * gz),
                    s[:, 1] * (1 + s[:, 0] * gx) * (1 + s[:, 2] * gz),
                    s[:, 2] * (1 + s[:, 0] * gx) * (1 + s[:, 1] * gy),
                ], axis=1) / 8 * 2  # d/dx = 2 d/dxi for a unit cube
                B = np.zeros((6, 24))
                B[0, 0::3] = dN[:, 0]
                B[1, 1::3] = dN[:, 1]
                B[2, 2::3] = dN[:, 2]
                B[3, 0::3] = dN[:, 1]; B[3, 1::3] = dN[:, 0]
                B[4, 1::3] = dN[:, 2]; B[4, 2::3] = dN[:, 1]
                B[5, 0::3] = dN[:, 2]; B[5, 2::3] = dN[:, 0]
                KE += B.T @ D @ B * (0.5 ** 3)
    return KE, corners


def main():
    nx, ny, nz, iters = (int(a) for a in sys.argv[1:5])
    volfrac = float(sys.argv[5]) if len(sys.argv) > 5 else 0.3
    rmin = float(sys.argv[6]) if len(sys.argv) > 6 else 1.5
    penal, emin, move = 3.0, 1e-9, 0.2
    KE, corners = hex8_ke(0.3)

    node = lambda i, j, k: i + (nx + 1) * (j + (ny + 1) * k)
    nel, ndof = nx * ny * nz, 3 * (nx + 1) * (ny + 1) * (nz + 1)
    edof = np.zeros((nel, 24), dtype=np.int64)
    centers = np.zeros((nel, 3))
    e = 0
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                ns = [node(i + c[0], j + c[1], k + c[2]) for c in corners]
                edof[e] = np.array([[3 * n, 3 * n + 1, 3 * n + 2] for n in ns]).ravel()
                centers[e] = (i + 0.5, j + 0.5, k + 0.5)
                e += 1
    iK = np.repeat(edof, 24, axis=1).ravel()
    jK = np.tile(edof, (1, 24)).ravel()

    fixed = np.array([3 * node(0, j, k) + d for k in range(nz + 1) for j in range(ny + 1) for d in range(3)])
    free = np.setdiff1d(np.arange(ndof), fixed)
    F = np.zeros(ndof)
    load_nodes = [node(nx, j, 0) for j in range(ny + 1)]
    for n in load_nodes:
        F[3 * n + 2] = -1.0 / len(load_nodes)

    # Hat filter, dense neighbourhood search (fine for reference sizes).
    rows, cols, vals = [], [], []
    r = int(np.ceil(rmin))
    idx = np.arange(nel).reshape(nz, ny, nx)
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                for dk in range(-r, r + 1):
                    for dj in range(-r, r + 1):
                        for di in range(-r, r + 1):
                            a, b, c = i + di, j + dj, k + dk
                            if 0 <= a < nx and 0 <= b < ny and 0 <= c < nz:
                                w = rmin - np.sqrt(di * di + dj * dj + dk * dk)
                                if w > 0:
                                    rows.append(idx[k, j, i]); cols.append(idx[c, b, a]); vals.append(w)
    H = coo_matrix((vals, (rows, cols)), shape=(nel, nel)).tocsr()
    Hs = np.asarray(H.sum(axis=1)).ravel()

    x = np.full(nel, volfrac)
    for it in range(1, iters + 1):
        xp = (H @ x) / Hs
        E = emin + xp ** penal * (1 - emin)
        K = coo_matrix(((KE.ravel()[None, :] * E[:, None]).ravel(), (iK, jK)), shape=(ndof, ndof)).tocsc()
        U = np.zeros(ndof)
        U[free] = spsolve(K[free][:, free], F[free])
        ue = U[edof]
        ce = np.einsum("ij,jk,ik->i", ue, KE, ue)
        c = float(F @ U)
        dc = -penal * xp ** (penal - 1) * (1 - emin) * ce
        dc = H.T @ (dc / Hs)
        dv = H.T @ (np.ones(nel) / Hs)
        lo, hi = 1e-30, 1e30
        while hi / lo > 1 + 1e-6:
            mid = np.sqrt(lo * hi)
            xn = np.clip(x * np.sqrt(np.maximum(-dc, 0) / (mid * np.maximum(dv, 1e-12))), np.maximum(x - move, 0), np.minimum(x + move, 1))
            if ((H @ xn) / Hs).mean() > volfrac:
                lo = mid
            else:
                hi = mid
        xn = np.clip(x * np.sqrt(np.maximum(-dc, 0) / (hi * np.maximum(dv, 1e-12))), np.maximum(x - move, 0), np.minimum(x + move, 1))
        change = float(np.abs(xn - x).max())
        x = xn
        print(json.dumps({"iter": it, "c": c, "vol": float(((H @ x) / Hs).mean()), "change": change}), flush=True)


if __name__ == "__main__":
    main()
