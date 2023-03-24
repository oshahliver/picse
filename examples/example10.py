from picse.interiors import planet_workbench
import pandas as pd

workbench = planet_workbench.Toolkit()

specs = {
    "ranges": {"mass": [0.01, 0.1], "x_FeS": [1e-6, 1e-6]},
    "sampling": {"mass": "log"},
    "core_segregation": "on",
}
inputs = workbench.sample_inputs(specs=specs, n_planets=1000, iteration_limit=1000)
df = pd.DataFrame(inputs)
print(df.head(10))
