import time
from picse.interiors import planet_workbench, planet_iterator, planet_creator
from picse.physicalparams import m_earth
from picse.materials import material

# Create an instance of the workbench toolkit
# Note. The EoS tables are loaded internally by the workbench
workbench = planet_workbench.Toolkit()

# Create an instance of the iterator toolkit
iterator = planet_iterator.Toolkit()

# Define core segregation pressure in Pa
total_mass = 0.1

t0 = time.time()
valid_compositions, pres_core_seg = workbench.create_core_mantle_composition(total_mass, num_models=100)
t1 = time.time()

# Compute outer core material fractions
outer_core_mat_fracs = [material.at2mat([1. - sum(vc[2:]), vc[2], vc[3], vc[4]]) for vc in valid_compositions]

for ocmf in outer_core_mat_fracs:
    print (ocmf, sum(ocmf))

# Calculate elapsed time
elapsed_time = t1 - t0

# Convert to minutes, seconds, and milliseconds
minutes = int(elapsed_time // 60)
seconds = int(elapsed_time % 60)
milliseconds = int((elapsed_time - minutes * 60 - seconds) * 1000)

# Output the valid compositions
# for i, comp in enumerate(valid_compositions):
#     print(f"Model {i+1}: P_CS = {1e-9*pres_core_seg[i]:.3g} GPa, Fe# = {comp[0]}, Si# = {comp[1]}, Core Composition: FeS = {comp[2]}, FeSi = {comp[3]}, FeO = {comp[4]}")

print(f"Elapsed Time: {minutes} min, {seconds} sec, {milliseconds} ms")

planetary_params = [{"M_surface_should":1.0,
                    "T_surface_should": 1800.0,
                    "Mg_number_should": 0.5162,
                    "Fe_number_mantle": vc[0],
                    "Si_number_mantle": vc[1],
                    "temperature_jumps":[0, 1800, 0, 0], # Temperature jumps across each layer transition
                    "contents": [[2], [2, 8, 10, 9], [4, 5], [6, 7]],
                    "fractions": [[1], ocmf, [.5, .5], [.5, .5]],
                    } for vc, ocmf in zip(valid_compositions, outer_core_mat_fracs)]

iterator_specs = [{
    "what": ["M_surface", "T_surface"],  # --> target properties
    "how": ["P_center", "T_center"],  # --> adjustable properties
    "val_should": [
        pp["M_surface_should"] * m_earth,
        pp["T_surface_should"],
    ],  # --> target values
    "predictor": ["linear", "linear"],  # --> no effect at this point
    "all_val_should_weights": [
        "log",
        "log",
    ],  # --> log or lin extrapolation for targets
    "all_howval_weights": ["exp", "exp"],  # --> exp or lin prediction for adjustables
    "acc": [1e-4, 1e-3],  # --> desired relative accuracies
    "iterationLimit": 20,  # --> max. number of iterations
    "deltaType": 0,  # --> mode for initial adjustment of adjustables
    "unpredictable": False,  # --> no effect at this point
} for pp in planetary_params]

run_params = {"layer_constraints": [1, 1, 3, 1], # Layer constraint type for each layer transition
              "core_segregation_type": 1,
              }

for pp in planetary_params:
    print (pp)

#######################################################################
# Model creation and execution
#######################################################################
planets = []
for pp, its in zip(planetary_params, iterator_specs):
    # Initialize a telluric planet instance with the specified properties
    pl = planet_creator.TelluricPlanet(planetary_params=pp, run_params = run_params)

    # Perform initial structure integration
    pl.construct()

    # Pass planet instance to iterator to match boundary conditions
    # NOTE. planetary objects passed to the iterator must be constructed!
    iterator.iterate(planet=pl, iterator_specs=its)
    planets.append(pl)

#######################################################################
# Model inspection
#######################################################################

# print fundamental planeatary properties to standard output
for pl in planets:
    pl.print(digits=4)