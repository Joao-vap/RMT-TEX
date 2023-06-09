import numpy as np

# ------------input-------------------------

n_particles = 20
nGroups = 5
steps = 100

# -------------functions-------------------

def ParticlesPerGroup(n_particles, nGroups):  
    if n_particles < nGroups:
        raise ValueError('n_particles must be greater than nGroups')
    particlesPerGroup = np.zeros(nGroups)
    for i in range(n_particles):
        particlesPerGroup[i%nGroups] += 1
    return [int(i) for i in particlesPerGroup]

def GUE(N, sigma = 1, mu = 0, escale = False):
    beta = 2
    A = np.random.randn(N,N) * sigma + 1j*np.random.randn(N,N)*sigma
    A = (A + A.T.conj())/2
    eig = np.linalg.eigvalsh(A)
    if escale:
        return (eig + mu)/(np.sqrt(N*beta))
    return eig + mu

# ------------parameters--------------------

particlesPerGroup = ParticlesPerGroup(n_particles, nGroups)
med = np.median(range(nGroups))
final_group_positions = [(i-med)*20 for i in range(nGroups)]

# ------------main--------------------------

eig_GUE = np.zeros(n_particles)

final_time = 1
dt = final_time/steps

Memory = np.zeros((steps+1,n_particles))
for step in range(steps+1):

    partial_time = step * dt
    eig_GUE = []

    sigma = 2*partial_time*(final_time-partial_time)/final_time
    for group in range(nGroups):
        mu = final_group_positions[group]*partial_time/final_time
        aux_eig_GUE = GUE(particlesPerGroup[group], sigma = sigma, mu = mu, escale = False)
        aux_ordered_eig_GUE = np.sort(aux_eig_GUE)

        if group == 0:
            eig_GUE = aux_ordered_eig_GUE
        else:
            while aux_ordered_eig_GUE[0] < eig_GUE[-1]:
                aux_eig_GUE = GUE(particlesPerGroup[group], sigma = sigma, mu = mu, escale = False)
                aux_ordered_eig_GUE = np.sort(aux_eig_GUE)
            
            eig_GUE = np.concatenate((eig_GUE, aux_ordered_eig_GUE))


    print(step)

    Memory[step,:] = eig_GUE

# ------------write Memory in file----------

np.savetxt('Memory.txt', Memory)

