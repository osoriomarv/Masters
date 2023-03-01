# # weight percent data of each silicate component
silicate_components = {'Si': 22.83, 'Mg': 25.20, 'Fe': 0.05, 'Al': 1.36, 'Ca': 2.98, 'Ni': 0.1639}

# weight percent data of each metal component
metal_components = {'Fe': 85.59, 'Ni': 5, 'O': 0, 'Si': 9.17}

# molar mass of each silicate component in g/mol
silicate_molar_mass = {'SiO2': 60.08, 'MgO': 40.31, 'FeO': 71.85, 'Al2O3': 101.96, 'CaO': 56.08, 'NiO': 74.69} 

# molar mass of each metal component in g/mol
Metal_molar_mass = {'Fe': 55.85, 'Ni': 58.69, 'O': 16, 'Si': 28.08}

# mass of the melt in g
mass = 1

# calculate the moles of each silicate component
silicate_moles = {}
for component, wt_percent in silicate_components.items():
    if component == 'Si':
        silicate_moles['SiO2'] = (wt_percent / 100) * mass / silicate_molar_mass['SiO2']
    elif component == 'Mg':
        silicate_moles['MgO'] = (wt_percent / 100) * mass / silicate_molar_mass['MgO']
    elif component == 'Fe':
        silicate_moles['FeO'] = (wt_percent / 100) * mass / silicate_molar_mass['FeO']
    elif component == 'Al':
        silicate_moles['Al2O3'] = (wt_percent / 100) * mass / silicate_molar_mass['Al2O3']
    elif component == 'Ca':
        silicate_moles['CaO'] = (wt_percent / 100) * mass / silicate_molar_mass['CaO']
    elif component == 'Ni':
        silicate_moles['NiO'] = (wt_percent / 100) * mass / silicate_molar_mass['NiO']

# calculate the moles of each metal component
metal_moles = {}
for component, wt_percent in metal_components.items():
    if component == 'Fe':
        metal_moles['Fe'] = (wt_percent / 100) * mass / Metal_molar_mass['Fe']
    elif component == 'Ni':
        metal_moles['Ni'] = (wt_percent / 100) * mass / Metal_molar_mass['Ni']
    elif component == 'O':
        metal_moles['O'] = (wt_percent / 100) * mass / Metal_molar_mass['O']
    elif component == 'Si':
        metal_moles['Si'] = (wt_percent / 100) * mass / Metal_molar_mass['Si']

# calculate the total moles of each component type
silicate_total_moles = sum(silicate_moles.values())

# calculate the total moles in metal
total_Metal_moles = sum(metal_moles.values())

# calculate the total moles
total_moles = silicate_total_moles + total_Metal_moles

# calculate the mol fraction of each silicate component
silicate_mol_fraction = {}
for component, m in silicate_moles.items():
    silicate_mol_fraction[component] = m / total_moles

# calculate the mol fraction of each metal component
Metal_mol_fraction = {}
for component, m in metal_moles.items():
    Metal_mol_fraction[component] = m / total_moles

print("Mol fraction of each silicate component:", silicate_mol_fraction)
print("Mol fraction of each metal component:", Metal_mol_fraction)


# # weight percent data of each silicate component
# silicate_components = {'SiO2': 44.51, 'MgO': 31.05, 'FeO': 17.52, 'Al2O3': 3.81, 'CaO': 3.1, 'NiO': 0}

# # weight percent data of each metal component
# metal_components = {'Fe': 91.04, 'Ni': 8.56, 'O': 0, 'Si': 0}
# # # weight percent data of each silicate component
# # silicate_components = {'SiO2': 48.84, 'MgO': 41.79, 'FeO': 0.06, 'Al2O3': 5.13, 'CaO': 4.17, 'NiO': 0}

# # # weight percent data of each metal component
# # metal_components = {'Fe': 85.59, 'Ni': 5, 'O': 0, 'Si': 9.17}

# # molar mass of each silicate component in g/mol
# silicate_molar_mass = {'SiO2': 60.08, 'MgO': 40.31, 'FeO': 71.85, 'Al2O3': 101.96, 'CaO': 56.08, 'NiO': 74.69}

# # molar mass of each metal component in g/mol
# metal_molar_mass = {'Fe': 55.85, 'Ni': 58.69, 'O': 16, 'Si': 28.08}

# # mass of the melt in g
# mass = 1

# # calculate the moles of each silicate component
# silicate_moles = {}
# for component, wt_percent in silicate_components.items():
#     silicate_moles[component] = (wt_percent / 100) * mass / silicate_molar_mass[component]

# # calculate the moles of each metal component
# metal_moles = {}
# for component, wt_percent in metal_components.items():
#     if component == 'O':
#         # skip the oxygen component, since it's not a metal
#         continue
#     metal_moles[component] = (wt_percent / 100) * mass / metal_molar_mass[component]

# # calculate the total moles of each component type
# silicate_total_moles = sum(silicate_moles.values())
# total_metal_moles = sum(metal_moles.values())

# # calculate the total moles
# total_moles = silicate_total_moles + total_metal_moles

# # calculate the mol fraction of each silicate component
# silicate_mol_fraction = {}
# for component, m in silicate_moles.items():
#     silicate_mol_fraction[component] = m / total_moles

# # calculate the mol fraction of each metal component
# metal_mol_fraction = {}
# for component, m in metal_moles.items():
#     metal_mol_fraction[component] = m / total_moles

# print("Mol fraction of each silicate component:", silicate_mol_fraction)
# print("Mol fraction of each metal component:", metal_mol_fraction)