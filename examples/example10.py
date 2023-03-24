from picse.interiors import planet_workbench
import pandas as pd

workbench = planet_workbench.Toolkit()

specs = {
    "ranges": {"mass": [0.5, 6]},
    "sampling": {"mass": "log", "xi_FeSi": "log", "xi_FeO": "log"},
    "core_segregation": "on",
}
inputs = workbench.sample_inputs(specs=specs, n_planets=10)
df = pd.DataFrame(inputs)
print(df.head(10))
