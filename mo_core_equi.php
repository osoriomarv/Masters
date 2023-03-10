<?php

require_once 'vendor/autoload.php';

use MathPHP\Functions\Special;

$P = 40;
$T_v = range(1000, 2500, round((2500-1000)*4));

// partitioning expressions
$KD_Si_p = [1.3, -13500, 0];
$KD_Ni_p = [0.46, 2700, -61];
$KD_O_p = [0.6, -3800, 22];

$KD_Si_v = Special::pow(10, $KD_Si_p[0] + $KD_Si_p[1] / $T_v + $KD_Si_p[2] * $P / $T_v);
$KD_Ni_v = Special::pow(10, $KD_Ni_p[0] + $KD_Ni_p[1] / $T_v + $KD_Ni_p[2] * $P / $T_v);
$KD_O_v = Special::pow(10, $KD_O_p[0] + $KD_O_p[1] / $T_v + $KD_O_p[2] * $P / $T_v);

// bulk composition,
// SiO2, MgO, FeO, Al2O3, CaO, NiO
$melt_wt = [48.84, 41.79, 0.06, 5.13, 4.17, 0.0];
$melt_mm = [28.19+16, 24.31+16, 55.85+16, 26.98+16*1.5, 40.08+16, 58.59+16];
$melt_mols_i = array_map(function($wt, $mm) {
    return $wt / $mm;
}, $melt_wt, $melt_mm);
$melt_mol_frac_i = array_map(function($mol) use ($melt_mols_i) {
    return $mol / array_sum($melt_mols_i);
}, $melt_mols_i);

// Fe, Ni, O, Si
$metal_wt = [85.89, 5.00, 0, 9.17];
$metal_mm = [55.85, 58.59, 16, 28.19];
$metal_mols_i = array_map(function($wt, $mm) {
    return $wt / $mm;
}, $metal_wt, $metal_mm);
$metal_mol_frac_i = array_map(function($mol) use ($metal_mols_i) {
    return $mol / array_sum($metal_mols_i);
}, $metal_mols_i);

$Fe_met_frac = 0.999;
$metal_silicate_ratio = $Fe_met_frac * $melt_mols_i[2] / ($metal_mols_i[0] - $Fe_met_frac * $metal_mols_i[0]);

$metal_mols_i_scaled = array_map(function($mol) use ($metal_silicate_ratio) {
    return $mol * $metal_silicate_ratio;
}, $metal_mols_i);

?>
