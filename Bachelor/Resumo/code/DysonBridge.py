import numpy as np

n_particles = 20
steps = 1000

def GUE(N, sigma = 1, mu = 0, escale = False):
    beta = 2
    A = np.random.randn(N,N) * sigma + 1j*np.random.randn(N,N)*sigma
    A = (A + A.T.conj())/2
    eig = np.linalg.eigvalsh(A)
    if escale:
        return (eig + mu)/(np.sqrt(N*beta))
    return eig + mu

eig_GUE = np.zeros(n_particles)

final_time = 1
dt = final_time/steps
Memory = np.zeros((steps+1,n_particles))

for step in range(steps+1):
    partial_time = step * dt
    sigma = 2*partial_time*(final_time-partial_time)/final_time
    mu = 0
    eig_GUE = GUE(n_particles, sigma = sigma, mu = mu, escale = False)
    eig_order = np.argsort(eig_GUE)
    Memory[step,:] = eig_GUE[eig_order]

np.savetxt('Memory.txt', Memory)
