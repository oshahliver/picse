from pics.materials import material

mix_specs = {"contents": [3, 5], "fractions": [0.25, 0.75], "temp": 300.0, "pres": 2e10}
mix = material.Mixture(specs=mix_specs)
mix.compute()
mix.update(temp=350.0, pres=3e10)
mix.print()
