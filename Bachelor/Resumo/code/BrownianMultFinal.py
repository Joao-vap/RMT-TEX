import numpy as np

n_particles = 100
n_groups = 2
steps = 100

def Mod_GUE(U, S, shift, N, sigma = 1, mu = 0, escale = False):
    beta = 2
    A = sigma * np.random.randn(N,N,) + sigma *1j*np.random.randn(N,N)
    A = (A + A.T.conj())/2
    A = A + (U @ S @ U.T.conj())*shift
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

def sigma_t(t):
    return 2*t*(1-t)

U, S = USU(n_groups, n_particles)
eig_GUE = np.zeros(n_particles)
dt = 1/steps

Memory = np.zeros((steps+1,n_particles))

for step in range(steps+1):
    print(step)
    partial_time = step * dt
    sigma = sigma_t(partial_time)
    shift = partial_time
    mu = 0
    eig_GUE = Mod_GUE(U, S, shift, n_particles, sigma = sigma, mu = mu, escale = False)
    eig_order = np.argsort(eig_GUE)
    Memory[step,:] = eig_GUE[eig_order]

# write Memory in file
np.savetxt('Memory.txt', Memory)
