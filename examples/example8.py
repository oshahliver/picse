from picse.materials import material
import time
import multiprocessing
from picse.interiors import planet_workbench

# Create an instance of the workbench toolkit
# Note. The EoS tables are loaded internally by the workbench
workbench = planet_workbench.Toolkit()

# Define core segregation pressure in Pa
total_mass = 0.1

# Number of models to generate
num_models = 10000  # Adjust this number as needed

t0 = time.time()
valid_compositions, pres_core_seg = workbench.create_core_mantle_composition(total_mass, num_models=num_models)
t1 = time.time()

# Calculate elapsed time
elapsed_time = t1 - t0

# Convert to minutes, seconds, and milliseconds
minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)
milliseconds = int((elapsed_time - minutes * 60 - seconds) * 1000)

# Output the valid compositions
for i, comp in enumerate(valid_compositions):
    print(f"Model {i+1}: P_CS = {1e-9*pres_core_seg[i]:.3g} GPa, Fe# = {comp[0]}, Si# = {comp[1]}, Core Composition: FeS = {comp[2]}, FeSi = {comp[3]}, FeO = {comp[4]}")

print(f"Elapsed Time: {minutes} min, {seconds} sec, {milliseconds} ms")
