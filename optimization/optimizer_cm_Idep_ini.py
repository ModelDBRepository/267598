from Genetic_optimization_cm_Idep_ini import runsimulation
from joblib import Parallel, delayed
import multiprocessing
import time
t = time.time()

prova=['neuron1']
num_cores = len(prova)

Ltt = Parallel(n_jobs=num_cores, verbose=50)(delayed(runsimulation)(prova[i],False) for i in range(len(prova)))

elapsed = time.time() - t
