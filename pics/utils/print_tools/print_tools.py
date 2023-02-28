from tabulate import tabulate
from pics.physicalparams import r_earth, m_earth, mH2O, material_list_fort
from pics.utils import functionTools as ftool


def print_planet(self, style=0, digits=3):
    """Prints out a simple overview of all relevant planetary parameters."""
    show_parameters = self.__dict__
    print("\n####################\nPlanet Label: {}\n".format(self.label))

    if style == 0:
        print("=======================================")
        print("Planet properties:")

        try:
            print("\n--------------------------------")
            print("Major parameters:")
            print("--------------------------------")
            print(
                "R_surface_is [R_earth]:",
                round(self.R_surface_is / r_earth, digits),
                "\nM_surface_is [M_earth]:",
                round(self.M_surface_is / m_earth, digits),
                "\nmean density [gcc]:",
                round(self.mean_density / 1000, digits),
                "\nT_surface_is [K]:",
                round(self.T_surface_is, digits),
                "\nT_surface_should [K]:",
                round(self.T_surface_should, digits),
                "\nT_center [K]:",
                round(self.T_center, digits),
                "\nP_surface_is [bar]:",
                round(self.P_surface_is * 1.0e-5, digits),
                "\nP_surface_should [bar]:",
                round(self.P_surface_should * 1.0e-5, digits),
                "\nP_center [GPa]:",
                round(self.P_center * 1.0e-9, digits),
                "\nMOI factor:",
                round(
                    self.moment_of_inertia_is
                    / (self.M_surface_is * self.R_surface_is**2),
                    digits,
                ),
                "\n\nMg_number_should:",
                round(self.Mg_number_should, digits),
                "\nMg_number_is:",
                round(self.Mg_number_is, digits),
                # "\nSi_number_should:",
                # round(self.Si_number_should, digits),
                "\nSi_number_is:",
                round(self.Si_number_is, digits),
            )

            try:
                print(
                    "M_H2O [wt%]:",
                    round(
                        (self.H2O_count + self.H_count / 2.0)
                        * mH2O
                        / self.M_surface_is
                        * 100,
                        digits,
                    ),
                )
                print(
                    "M_H2O core [wt%]:",
                    round(self.H_count / 2.0 * mH2O / self.M_surface_is * 100, digits),
                )
                print(
                    "M_H2O mantle [wt%]:",
                    round(self.H2O_count * mH2O / self.M_surface_is * 100, digits),
                )
                print("Ocean fraction is:", round(10**self.ocean_fraction_is, digits))
                print(
                    "Ocean fraction should:",
                    round(10**self.ocean_fraction_should, digits),
                )
                print(
                    "Core mass fraction:",
                    round(self.M_core_is / self.M_surface_is * m_earth, digits),
                )
                print(
                    "Core radius [km]:",
                    round(self.layer_properties[1]["R_outer"] / 1000, digits),
                )
            except ZeroDivisionError:
                print("M_H2O [wt%]: NaN")

        except TypeError:
            print("WARNING: Type Error in Planet.prt()")

        print("\n--------------------------------")
        print("Layer overview:")
        print("--------------------------------")

        dat = []
        for i in range(len(self.layer_properties)):
            lay = self.layer_properties[i]
            material_str = ""
            for c in range(len(self.contents[i])):
                frac = str(round(self.fractions[i][c] * 100, 1)) + "% "
                material_str += frac + material_list_fort[self.contents[i][c] - 1]
                if c < len(self.contents[i]) - 1:
                    material_str += ", "

            dat.append(
                [
                    i,
                    material_str,
                    ftool.scinot(lay["R_outer"] / r_earth, digits=digits),
                    ftool.scinot(lay["indigenous_mass"] / m_earth, digits=digits),
                    ftool.scinot(lay["P_outer"] * 1.0e-9, digits=digits),
                    ftool.scinot(lay["T_outer"], digits=digits),
                    ftool.scinot(lay["rho_outer"], digits=digits),
                ]
            )

        tabl = tabulate(
            dat,
            headers=[
                "Layer",
                "Contents",
                "R [R_e]",
                "m [M_e]",
                "P [GPa]",
                "T [K]",
                "rho [kg m-3]",
            ],
        )

        print()
        print(f"{tabl}")
        print()
