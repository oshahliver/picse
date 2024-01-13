import numpy as np
import random
from picse.materials import material
import time

t0 = time.time()

def sample_log_uniform(low, high):
    """Sample a value uniformly on a logarithmic scale."""
    return np.exp(random.uniform(np.log(low), np.log(high)))

# Define core segregation pressure in Pa
pres_core_seg = 5e10

# Number of models to generate
num_models = 10000  # Adjust this number as needed

# =====================================================================
# Array to store valid core-mantle composition pairs
# valid_compositions = []

# for _ in range(num_models):
#     valid_composition = False
#     while not valid_composition:
#         # Randomly sample core composition components on a log scale
#         xiFeS = sample_log_uniform(0.001, 0.5)  # Avoid zero in log scale
#         xiFeSi = sample_log_uniform(0.001, 0.3)
#         xiFeO = sample_log_uniform(0.001, 0.2)

#         xi_all_core = [0.0, xiFeS, xiFeSi, xiFeO]  # xiFe will be computed automatically
#         at = material.mat2at_core(xi=xi_all_core, xiH=0)

#         # Compute liquidus of pyrolite and mantle composition
#         temp_core_seg = material.temp_liquidus_pyrolite(pres_core_seg)
#         Fe_number_mantle, Si_number_mantle = material.calc_equilibrium_mantle_comp(pres_core_seg, xi_all_core)

#         # Check if the composition is within the feasible range
#         simmin = material.Si_number_mantle_min(Fe_number_mantle)
#         simmax = material.Si_number_mantle_max(Fe_number_mantle)
#         femmin = 0.0
#         femmax = 0.5

#         if simmin <= Si_number_mantle <= simmax and femmin <= Fe_number_mantle <= femmax:
#             valid_composition = True
#             valid_compositions.append((Fe_number_mantle, Si_number_mantle, xiFeS, xiFeSi, xiFeO))

# t1 = time.time()

# # Calculate elapsed time
# elapsed_time = t1 - t0

# # Convert to minutes, seconds, and milliseconds
# minutes = int(elapsed_time // 60)
# seconds = int(elapsed_time % 60)
# milliseconds = int((elapsed_time - minutes * 60 - seconds) * 1000)

# print(f"Elapsed Time: {minutes} min, {seconds} sec, {milliseconds} ms")

# Output the valid compositions
# for i, comp in enumerate(valid_compositions):
#     print(f"Model {i+1}: Fe# = {comp[0]}, Si# = {comp[1]}, Core Composition: FeS = {comp[2]}, FeSi = {comp[3]}, FeO = {comp[4]}")
# =====================================================================

# import numpy as np
# from picse.materials import material

# t0 = time.time()

# def vectorized_log_uniform_sampling(low, high, size):
#     """Vectorized sampling of values uniformly on a logarithmic scale."""
#     return np.exp(np.random.uniform(np.log(low), np.log(high), size))

# # Initialize arrays to store core composition components and results
# xiFeS_samples = np.zeros(num_models)
# xiFeSi_samples = np.zeros(num_models)
# xiFeO_samples = np.zeros(num_models)
# Fe_number_mantle_samples = np.zeros(num_models)
# Si_number_mantle_samples = np.zeros(num_models)

# # Initialize an array to track which models are valid
# valid_models = np.zeros(num_models, dtype=bool)

# while not np.all(valid_models):
#     # Indices of models that need resampling
#     invalid_indices = np.where(~valid_models)[0]

#     # Resample invalid models
#     xiFeS_samples[invalid_indices] = vectorized_log_uniform_sampling(0.001, 0.5, len(invalid_indices))
#     xiFeSi_samples[invalid_indices] = vectorized_log_uniform_sampling(0.001, 0.3, len(invalid_indices))
#     xiFeO_samples[invalid_indices] = vectorized_log_uniform_sampling(0.001, 0.2, len(invalid_indices))

#     # Compute mantle compositions for resampled models
#     for i in invalid_indices:
#         xi_all_core = [0.0, xiFeS_samples[i], xiFeSi_samples[i], xiFeO_samples[i]]
#         at = material.mat2at_core(xi=xi_all_core, xiH=0)
#         temp_core_seg = material.temp_liquidus_pyrolite(pres_core_seg)
#         Fe_number_mantle, Si_number_mantle = material.calc_equilibrium_mantle_comp(pres_core_seg, xi_all_core)
#         Fe_number_mantle_samples[i] = Fe_number_mantle
#         Si_number_mantle_samples[i] = Si_number_mantle

#         # Check if the composition is valid
#         if (material.Si_number_mantle_min(Fe_number_mantle) <= Si_number_mantle <= material.Si_number_mantle_max(Fe_number_mantle)) and (0.0 <= Fe_number_mantle <= 0.5):
#             valid_models[i] = True

# t1 = time.time()

# # Calculate elapsed time
# elapsed_time = t1 - t0

# # Convert to minutes, seconds, and milliseconds
# minutes = int(elapsed_time // 60)
# seconds = int(elapsed_time % 60)
# milliseconds = int((elapsed_time - minutes * 60 - seconds) * 1000)

# print(f"Elapsed Time: {minutes} min, {seconds} sec, {milliseconds} ms")

# # Output the valid compositions
# # valid_compositions = np.column_stack((Fe_number_mantle_samples, Si_number_mantle_samples, xiFeS_samples, xiFeSi_samples, xiFeO_samples))
# # for i, comp in enumerate(valid_compositions):
# #     print(f"Model {i+1}: Fe# = {comp[0]}, Si# = {comp[1]}, Core Composition: FeS = {comp[2]}, FeSi = {comp[3]}, FeO = {comp[4]}")
# =====================================================================
import numpy as np
import random
import multiprocessing
from picse.materials import material

t0 = time.time()

def compute_model_wrapper(_):
    """Wrapper function for compute_model to be used with multiprocessing.
    
    Args:
        dummy_arg: A placeholder argument, required for compatibility with multiprocessing's map function.
    """
    return material.compute_model(pres_core_seg)

# Use multiprocessing pool
with multiprocessing.Pool(processes=multiprocessing.cpu_count()) as pool:
    valid_compositions = pool.map(compute_model_wrapper, range(num_models))

pool.close()
pool.join()

t1 = time.time()

# Calculate elapsed time
elapsed_time = t1 - t0

# Convert to minutes, seconds, and milliseconds
minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)
milliseconds = int((elapsed_time - minutes * 60 - seconds) * 1000)

# Output the valid compositions
for i, comp in enumerate(valid_compositions):
    print(f"Model {i+1}: Fe# = {comp[0]}, Si# = {comp[1]}, Core Composition: FeS = {comp[2]}, FeSi = {comp[3]}, FeO = {comp[4]}")

print(f"Elapsed Time: {minutes} min, {seconds} sec, {milliseconds} ms")
