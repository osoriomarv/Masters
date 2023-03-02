def element_distribution(Si, Ca, Al, Mg, O, Fe, Ni):
    # Calculate the total weight in mg/g
    total_weight = Si + Ca + Al + Mg + O + Fe + Ni

    # Calculate the amount of oxygen needed for each oxide in mg/g
    MgO = Mg
    Al2O3 = Al * 1.5
    CaO = Ca
    remaining_O = O - (MgO + Al2O3 + CaO)

    # Calculate the amount of SiO2 that can be formed from the remaining oxygen
    SiO2 = min(remaining_O, Si)
    remaining_O -= SiO2
    Si -= SiO2

    # Calculate the amount of Si and free oxygen
    free_O = remaining_O
    remaining_O -= free_O
    free_Si = max(Si - SiO2, 0)

    # Calculate the weight percentage of each oxide
    wt_MgO = 100 * MgO / total_weight
    wt_Al2O3 = 100 * Al2O3 / total_weight
    wt_CaO = 100 * CaO / total_weight
    wt_SiO2 = 100 * SiO2 / total_weight
    wt_free_O = 100 * free_O / total_weight

    # Calculate the weight percentage of each element
    wt_Si = 100 * Si / total_weight
    wt_Ca = 100 * Ca / total_weight
    wt_Al = 100 * Al / total_weight
    wt_Mg = 100 * Mg / total_weight
    wt_Fe = 100 * Fe / total_weight
    wt_Ni = 100 * Ni / total_weight

    # Calculate the molar fraction of each element
    molar_masses = {'Si': 28.085, 'Ca': 40.078, 'Al': 26.982, 'Mg': 24.305,
                    'O': 15.999, 'Fe': 55.845, 'Ni': 58.693}
    total_moles = sum([elem_weight / molar_masses[elem] for elem, elem_weight in
                       {'Si': Si, 'Ca': Ca, 'Al': Al, 'Mg': Mg, 'O': MgO + Al2O3 + CaO + SiO2 + free_O,
                        'Fe': Fe, 'Ni': Ni}.items()])
    mol_Si = (Si / molar_masses['Si']) / total_moles
    mol_Ca = (Ca / molar_masses['Ca']) / total_moles
    mol_Al = (Al / molar_masses['Al']) / total_moles
    mol_Mg = (Mg / molar_masses['Mg']) / total_moles
    mol_O = ((MgO + Al2O3 + CaO + SiO2 + free_O) / molar_masses['O']) / total_moles
    mol_Fe = (Fe / molar_masses['Fe']) / total_moles
    mol_Ni = (Ni / molar_masses['Ni']) / total_moles
    mol_MgO = (MgO / molar_masses['Mg'] + 0.5 * mol_O) / total_moles
    mol_Al2O3 = (Al2O3 / molar_masses['Al'] + 1.5*mol_O) / total_moles
    mol_CaO = (CaO / molar_masses['Ca'] + mol_O) / total_moles 
    mol_SiO2 = (SiO2 / molar_masses['Si'] + 2 * mol_O) / total_moles

    # Return the results as a dictionary
    return {'MgO': {'wt%': wt_MgO, 'mol_frac': mol_MgO, 'mg/g': MgO},
            'Al2O3': {'wt%': wt_Al2O3, 'mol_frac': mol_Al2O3, 'mg/g': Al2O3},
            'CaO': {'wt%': wt_CaO, 'mol_frac': mol_CaO, 'mg/g': CaO},
            'SiO2': {'wt%': wt_SiO2, 'mol_frac': mol_SiO2, 'mg/g': SiO2},
            'free_O': {'wt%': wt_free_O, 'mg/g': free_O},
            'Si': {'wt%': wt_Si, 'mol_frac': mol_Si, 'mg/g': Si},
            'Ca': {'wt%': wt_Ca, 'mol_frac': mol_Ca, 'mg/g': Ca},
            'Al': {'wt%': wt_Al, 'mol_frac': mol_Al, 'mg/g': Al},
            'Mg': {'wt%': wt_Mg, 'mol_frac': mol_Mg, 'mg/g': Mg},
            'Fe': {'wt%': wt_Fe, 'mol_frac': mol_Fe, 'mg/g': Fe},
            'Ni': {'wt%': wt_Ni, 'mol_frac': mol_Ni, 'mg/g': Ni},
            'free_Si': {'mg/g': free_Si}}

if __name__ == "__main__": 
    result = element_distribution(186, 10.1, 10.5, 141, 310, 220, 13)
    for elem, data in result.items():
        print(elem + ':')
        print('\tmg/g:', data['mg/g'])
        if 'mol_frac' in data:
            print('\tmol fraction:', data['mol_frac'])
        print('\twt%:', data['wt%'])