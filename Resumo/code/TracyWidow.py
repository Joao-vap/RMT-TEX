import numpy as np

n_particles = 1000
n_groups = 2

def Mod_GUE(S, shift, N, sigma = 1, mu = 0, escale = False):
    beta = 2
    A = sigma * np.random.randn(N,N,) + sigma *1j*np.random.randn(N,N)
    A = (A + A.T.conj())/2
    A = A + S*shift
    eig = np.linalg.eigvalsh(A)
    if escale:
        return (eig*sigma + mu)/(np.sqrt(N*beta))
    return eig

def USU(nGroups, n_particles):
    # U is the unitary matrix with the e1, e2, ..., eN vectors as columns
    # S is the diagonal matrix with the eigenvalues

    U = np.zeros((n_particles, n_particles), dtype = complex)
    S = np.zeros((n_particles, n_particles), dtype = complex)

    # fill U diagonal with 1's
    for i in range(n_particles):
        U[i,i] = 1
    
    # fill S diagonal with nGroups different values centered at 0
    divider = n_particles/nGroups
    shift = (nGroups-1)/2
    for i in range(n_particles):
        S[i,i] = (int(i/divider) - shift) * 10

    return U, S

def isHermitian(A):
    return np.allclose(A, A.T.conj())

def sigma_t(t):
    return 2*t*(1-t)

U, S = USU(n_groups, n_particles)
eig_GUE = np.zeros(n_particles)

exps = 1000

Memory = np.zeros((exps, n_particles))

for exp in range(exps):
    partial_time = 0.5
    sigma = sigma_t(partial_time)
    shift = partial_time
    mu = 0
    eig_GUE = Mod_GUE(S, shift, n_particles, sigma = sigma, mu = mu, escale = False)
    eig_order = np.argsort(eig_GUE)
    Memory[exp,:] = eig_GUE[eig_order]

    print(exp)

# write Memory in file
np.savetxt('Memory_test.txt', Memory)
