import numpy as np
import matplotlib.pyplot as plt

P = 10
T_v = np.linspace(2000, 2000, 1)

# partitioning expressions
KD_Si_p = [1.3, -13500, 0]
KD_Ni_p = [0.46, 2700, -61]
KD_O_p = [0.6, -3800, 22]

# Solve for the Fischer values - vectorized
KD_Si_v = 10 ** (KD_Si_p[0] + KD_Si_p[1] / T_v + KD_Si_p[2] * P / T_v)
KD_Ni_v = 10 ** (KD_Ni_p[0] + KD_Ni_p[1] / T_v + KD_Ni_p[2] * P / T_v)
KD_O_v = 10 ** (KD_O_p[0] + KD_O_p[1] / T_v + KD_O_p[2] * P / T_v)

# bulk composition
# SiO2, MgO, FeO, Al2O3, CaO, NiO
melt_wt = np.array([48.84, 41.79, 0.06, 5.13, 4.17, 0.0])
melt_mm = np.array([28.19 + 16, 24.31 + 16, 55.85 + 16, 26.98 + 16 * 1.5, 40.08 + 16, 58.69 + 16])
melt_mols_i = melt_wt / melt_mm
melt_mol_frac_i = melt_mols_i / np.sum(melt_mols_i)

# Fe, Ni, O, Si
metal_wt = np.array([85.89, 5.00, 0, 9.17])
metal_mm = np.array([55.85, 58.69, 16, 28.19])
metal_mols_i = metal_wt / metal_mm
metal_mol_frac_i = metal_mols_i / np.sum(metal_mols_i)

Fe_met_frac = 0.999
metal_silicate_ratio = Fe_met_frac * melt_mols_i[2] / (metal_mols_i[0] - Fe_met_frac * metal_mols_i[0])

metal_mols_i_scaled = metal_mols_i * metal_silicate_ratio

# bulk cation, Si, Mg, Fe, Al, Ca, Ni
bulk_cat = np.array([melt_mols_i[0] + metal_mols_i_scaled[3], melt_mols_i[1],
                     melt_mols_i[2] + metal_mols_i_scaled[0], melt_mols_i[3],
                     melt_mols_i[4], metal_mols_i_scaled[1]])

# bulk O
bulk_O = melt_mols_i[0] * 2 + melt_mols_i[1] + melt_mols_i[2] + melt_mols_i[3] * 1.5 + melt_mols_i[4]

bulk_comp = np.concatenate((bulk_cat, [bulk_O]))

FeO_guess_v = np.linspace(0.001, 0.1, int(1e3))

# x = FeO; y = NiO, z = SiO2, u = MgO, m = Al1.5O, n = CaO
# a = Fe, b = Ni, c = O, d
FeO_i = melt_mols_i[2]
NiO_i = melt_mols_i[5]
SiO2_i = melt_mols_i[0]
MgO_i = melt_mols_i[1]
AlO1_5_i = melt_mols_i[3]
CaO_i = melt_mols_i[4]

Fe_i = metal_mols_i_scaled[0]
Ni_i = metal_mols_i_scaled[1]
O_i = metal_mols_i_scaled[2]
Si_i = metal_mols_i_scaled[3]

FeO_mol_frac = np.zeros(len(FeO_guess_v))
NiO_mol_frac = np.zeros(len(FeO_guess_v))
SiO2_mol_frac = np.zeros(len(FeO_guess_v))
Fe_mol_frac = np.zeros(len(FeO_guess_v))
Ni_mol_frac = np.zeros(len(FeO_guess_v))
O_mol_frac = np.zeros(len(FeO_guess_v))
Si_mol_frac = np.zeros(len(FeO_guess_v))
KD_O_guess = np.zeros(len(FeO_guess_v))
iternum = np.arange(len(FeO_guess_v))
min_pos_2_logger = np.zeros((len(T_v), len(FeO_guess_v) - 1))
min_pos_logger = np.zeros(len(T_v), dtype=int)
FeO_mol_frac_pick = np.zeros(len(T_v))
NiO_mol_frac_pick = np.zeros(len(T_v))
SiO2_mol_frac_pick = np.zeros(len(T_v))
Fe_mol_frac_pick = np.zeros(len(T_v))
Ni_mol_frac_pick = np.zeros(len(T_v))
O_mol_frac_pick = np.zeros(len(T_v))
Si_mol_frac_pick = np.zeros(len(T_v))
KD_O_frac_er = np.zeros(len(T_v))
melt_mols_logger = np.zeros(len(T_v))
metal_mols_logger = np.zeros(len(T_v))
KD_O_guess_logger = np.zeros((len(FeO_guess_v), len(T_v)))

for q in range(len(T_v)):
    T = T_v[q]
    for w in range(len(FeO_guess_v)):
        FeO_guess = FeO_guess_v[w]
        NiO = FeO_guess * (NiO_i + Ni_i) / ((FeO_i + Fe_i - FeO_guess) * KD_Ni_v[q] + FeO_guess)
        Fe = FeO_i + Fe_i - FeO_guess
        Ni = NiO_i + Ni_i - NiO
        alpha = SiO2_i + Si_i
        gamma = Fe + Ni + FeO_i + NiO_i + 3 * SiO2_i + O_i - FeO_guess - NiO + Si_i
        sigma = FeO_guess + NiO + MgO_i + AlO1_5_i + CaO_i

        quad_a = (3 * FeO_guess ** 2 - Fe ** 2 * KD_Si_v[q])
        quad_b = -(gamma * FeO_guess ** 2 + 3 * alpha * FeO_guess ** 2 + Fe ** 2 * sigma * KD_Si_v[q])
        quad_c = alpha * gamma * FeO_guess ** 2
        root_solve = np.roots([quad_a, quad_b, quad_c])
        SiO2 = root_solve[1]
        O = FeO_i + NiO_i + SiO2_i * 2 + O_i + FeO_guess - NiO - 2 * SiO2
        Si = SiO2_i + Si_i - SiO2

        # x = FeO; y = NiO, z = SiO2
        # u = MgO, m = Al1.5O, n = CaO
        # a = Fe, b = Ni, c = O, d = Si
        melt_mols = FeO_guess + NiO + SiO2 + MgO_i + AlO1_5_i + CaO_i
        metal_mols = Fe + Ni + O + Si

        FeO_mol_frac[w] = FeO_guess / melt_mols
        NiO_mol_frac[w] = NiO / melt_mols
        SiO2_mol_frac[w] = SiO2 / melt_mols
        Fe_mol_frac[w] = Fe / metal_mols
        Ni_mol_frac[w] = Ni / metal_mols
        O_mol_frac[w] = O / metal_mols
        Si_mol_frac[w] = Si / metal_mols

        melt_mols_logger[q] = melt_mols
        metal_mols_logger[q] = metal_mols

        KD_O_guess[w] = Fe_mol_frac[w] * O_mol_frac[w] / FeO_mol_frac[w]
        iternum[w] = w
        KD_O_guess_logger[w, q] = Fe_mol_frac[w] * O_mol_frac[w] / FeO_mol_frac[w]

e = 1
for w in range(len(FeO_guess_v) - 1):
    if abs(KD_O_guess[w] - KD_O_v[q]) > abs(KD_O_guess[w + 1] - KD_O_v[q]):
        min_pos_2_logger[q, e] = w + 1
        min_pos_2 = w
    else:
        e += 1

min_pos = np.argmin(abs(KD_O_guess - KD_O_v[q]))
min_pos_logger[q] = min_pos

FeO_mol_frac_pick[q] = FeO_mol_frac[min_pos]
NiO_mol_frac_pick[q] = NiO_mol_frac[min_pos]
SiO2_mol_frac_pick[q] = SiO2_mol_frac[min_pos]
Fe_mol_frac_pick[q] = Fe_mol_frac[min_pos]
Ni_mol_frac_pick[q] = Ni_mol_frac[min_pos]
O_mol_frac_pick[q] = O_mol_frac[min_pos]
Si_mol_frac_pick[q] = Si_mol_frac[min_pos]
KD_O_frac_er[q] = (KD_O_guess[min_pos] - KD_O_v[q]) / KD_O_v[q]
FeO_guess_v = np.linspace(FeO_guess_v[min_pos] * 0.98, FeO_guess_v[min_pos] * 1.04, int(1e3))

fig1 = plt.figure(1)
plt.plot(FeO_guess_v, np.log10(KD_O_guess), ':k', [min(FeO_guess_v), max(FeO_guess_v)], [np.log10(KD_O_v[q]), np.log10(KD_O_v[q])], '-k')
plt.xlabel('FeO guess')
plt.ylabel('log10(KDO)')
plt.legend(['guess', 'Fischer'])
plt.show()


# figure 2
fig2 = plt.figure(2)
plt.subplot(3, 3, 1)
plt.plot(FeO_guess_v, SiO2_mol_frac, '-k', max(FeO_guess_v), melt_mol_frac_i[0], 'ok')
plt.title('SiO2')
plt.subplot(3, 3, 2)
plt.plot(FeO_guess_v, FeO_mol_frac, '-k', max(FeO_guess_v), melt_mol_frac_i[2], 'ok')
plt.title('FeO')
plt.subplot(3, 3, 3)
plt.plot(FeO_guess_v, NiO_mol_frac, '-k', max(FeO_guess_v), melt_mol_frac_i[5], 'ok')
plt.title('NiO')
plt.subplot(3, 3, 4)
plt.plot(FeO_guess_v, Fe_mol_frac, '-k', max(FeO_guess_v), metal_mol_frac_i[0], 'ok')
plt.title('Fe')
plt.subplot(3, 3, 5)
plt.plot(FeO_guess_v, Ni_mol_frac, '-k', max(FeO_guess_v), metal_mol_frac_i[1], 'ok')
plt.title('Ni')
plt.subplot(3, 3, 6)
plt.plot(FeO_guess_v, O_mol_frac, '-k', max(FeO_guess_v), metal_mol_frac_i[2], 'ok')
plt.title('O')
plt.subplot(3, 3, 7)
plt.plot(FeO_guess_v, Si_mol_frac, '-k', max(FeO_guess_v), metal_mol_frac_i[3], 'ok')
plt.title('Si')
plt.show()

# figure 3
fig3 = plt.figure(3)
plt.subplot(3, 3, 1)
plt.plot(T_v, SiO2_mol_frac_pick, '-k', [min(T_v), max(T_v)], [melt_mol_frac_i[0], melt_mol_frac_i[0]], ':k')
plt.title('SiO2')
plt.subplot(3, 3, 2)
plt.plot(T_v, FeO_mol_frac_pick, '-k', [min(T_v), max(T_v)], [melt_mol_frac_i[2], melt_mol_frac_i[2]], ':k')
plt.title('FeO')
plt.subplot(3, 3, 3)
plt.plot(T_v, NiO_mol_frac_pick, '-k', [min(T_v), max(T_v)], [melt_mol_frac_i[5], melt_mol_frac_i[5]], ':k')
plt.title('NiO')
plt.subplot(3, 3, 4)
plt.plot(T_v, Fe_mol_frac_pick, '-k', [min(T_v), max(T_v)], [metal_mol_frac_i[0], metal_mol_frac_i[0]], ':k')
plt.title('Fe')
plt.subplot(3, 3, 5)
plt.plot(T_v, Ni_mol_frac_pick, '-k', [min(T_v), max(T_v)], [metal_mol_frac_i[1], metal_mol_frac_i[1]], ':k')
plt.title('Ni')
plt.subplot(3, 3, 6)
plt.plot(T_v, O_mol_frac_pick, '-k', [min(T_v), max(T_v)], [metal_mol_frac_i[2], metal_mol_frac_i[2]], ':k')
plt.title('O')
plt.subplot(3, 3, 7)
plt.plot(T_v, Si_mol_frac_pick, '-k', [min(T_v), max(T_v)], [metal_mol_frac_i[3], metal_mol_frac_i[3]], ':k')
plt.title('Si')
plt.show()