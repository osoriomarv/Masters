def asteroid_composition(Si_wt, Ca_wt, Al_wt, Mg_wt, O_wt, Fe_wt, Ni_wt, Si_ppm, Ca_ppm, Al_ppm, Mg_ppm, Fe_ppm, Ni_ppm):
    # Calculate the total weight of the asteroid in grams
    total_weight = 100.0 / (Si_wt + Ca_wt + Al_wt + Mg_wt + O_wt + Fe_wt + Ni_wt) * 1000000.0

    # Convert the weight percentages to masses in grams
    Si = Si_wt / 100.0 * total_weight
    Ca = Ca_wt / 100.0 * total_weight
    Al = Al_wt / 100.0 * total_weight
    Mg = Mg_wt / 100.0 * total_weight
    O = O_wt / 100.0 * total_weight
    Fe = Fe_wt / 100.0 * total_weight
    Ni = Ni_wt / 100.0 * total_weight

    # Convert the ppm values to masses in grams
    Si += Si_ppm / 1000000.0 * total_weight
    Ca += Ca_ppm / 1000000.0 * total_weight
    Al += Al_ppm / 1000000.0 * total_weight
    Mg += Mg_ppm / 1000000.0 * total_weight
    Fe += Fe_ppm / 1000000.0 * total_weight
    Ni += Ni_ppm / 1000000.0 * total_weight

    # Calculate the total number of oxygen atoms
    total_O = O - (Mg + Al + Ca)

    # Calculate the number of oxygen atoms required to satisfy MgO
    MgO_O = Mg
    total_O -= MgO_O

    # Calculate the number of oxygen atoms required to satisfy Al2O3
    Al2O3_O = Al * 3
    total_O -= Al2O3_O

    # Calculate the number of oxygen atoms required to satisfy CaO
    CaO_O = Ca
    total_O -= CaO_O

    # Calculate the number of oxygen atoms available for SiO2
    SiO2_O = max(total_O, 0)

    # Calculate the number of Si atoms available
    Si_left = Si - SiO2_O

    # Report the composition
    print("MgO:", Mg)
    print("Al2O3:", Al)
    print("CaO:", Ca)
    print("SiO2:", SiO2_O)
    print("Si:", Si_left)
    print("Fe:", Fe)
    print("Ni:", Ni)
    print("Free O:", max(total_O, 0))

if __name__ == '__main__':
    # Enstatite Comp
    
    asteroid_composition(2.3, 4.5, 6.7, 8.9, 70.1, 0.01, 0.02, 123.4, 567.8, 910.1, 111.2, 33.4, 55.6)