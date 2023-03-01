#weight percent of Ni = (28 / 10,000) * 58.69 * 100 / 100 = 0.1639

# weight percent data of each oxide component
oxide_components = {'SiO2': 48.84, 'MgO': 41.79, 'FeO': 0.06, 'Al2O3': 5.13, 'CaO': 4.17}

# molar mass of each oxide component in g/mol
oxide_molar_mass = {'SiO2': 60.08, 'MgO': 40.31, 'FeO': 71.85, 'Al2O3': 101.96, 'CaO': 56.08}

# atomic mass of each element
atomic_mass = {'Si': 28.09, 'Mg': 24.31, 'Fe': 55.85, 'Al': 26.98, 'Ca': 40.08, 'Ni': 58.69}

# calculate the weight percent of each element
element_components = {}
for oxide, wt_percent in oxide_components.items():
    for element, mass in atomic_mass.items():
        if element in oxide:
            if element not in element_components:
                element_components[element] = 0
            element_components[element] += wt_percent * (mass / oxide_molar_mass[oxide])
element_components['Ni'] = 0.1639  # add the weight percent of Ni

print("Weight percent of each element:")
for element, wt_percent in element_components.items():
    print(element + ":", "{:.2f}".format(wt_percent) + " wt.%")
