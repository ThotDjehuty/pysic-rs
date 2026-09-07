//! Physical constants in SI units.
//!
//! All constants are provided via function accessors to enable future
//! support for unit system switching (SI / natural / CGS).

/// Speed of light in vacuum [m/s]
#[inline]
pub fn c() -> f64 { 299_792_458.0 }

/// Reduced Planck constant [J·s]
#[inline]
pub fn hbar() -> f64 { 1.054_571_817e-34 }

/// Planck constant [J·s]
#[inline]
pub fn h() -> f64 { 6.626_070_15e-34 }

/// Gravitational constant [m³·kg⁻¹·s⁻²]
#[inline]
pub fn G() -> f64 { 6.674_30e-11 }

/// Boltzmann constant [J/K]
#[inline]
pub fn k_B() -> f64 { 1.380_649e-23 }

/// Elementary charge [C]
#[inline]
pub fn e() -> f64 { 1.602_176_634e-19 }

/// Electron mass [kg]
#[inline]
pub fn m_e() -> f64 { 9.109_383_7015e-31 }

/// Proton mass [kg]
#[inline]
pub fn m_p() -> f64 { 1.672_621_92369e-27 }

/// Neutron mass [kg]
#[inline]
pub fn m_n() -> f64 { 1.674_927_49804e-27 }

/// Vacuum permittivity [F/m]
#[inline]
pub fn eps_0() -> f64 { 8.854_187_8128e-12 }

/// Vacuum permeability [H/m]
#[inline]
pub fn mu_0() -> f64 { 1.256_637_06212e-6 }

/// Fine-structure constant (dimensionless)
#[inline]
pub fn alpha() -> f64 { 7.297_352_5693e-3 }

/// Stefan-Boltzmann constant [W·m⁻²·K⁻⁴]
#[inline]
pub fn sigma_SB() -> f64 { 5.670_374_419e-8 }

/// Avogadro number [mol⁻¹]
#[inline]
pub fn N_A() -> f64 { 6.022_140_76e23 }

/// Rydberg constant [m⁻¹]
#[inline]
pub fn R_inf() -> f64 { 10_973_731.568_160 }

/// Compton wavelength of the electron [m]
#[inline]
pub fn lambda_C() -> f64 { 2.426_310_238_67e-12 }

/// Bohr radius [m]
#[inline]
pub fn a_0() -> f64 { 5.291_772_10903e-11 }

/// Classical electron radius [m]
#[inline]
pub fn r_e() -> f64 { 2.817_940_3262e-15 }

/// Hubble constant [km/s/Mpc] (approximate)
#[inline]
pub fn H_0() -> f64 { 67.4 }

/// Cosmological constant [m⁻²] (approximate, from Planck 2018)
#[inline]
pub fn Lambda() -> f64 { 1.105_6e-52 }

/// Planck length [m]
#[inline]
pub fn l_P() -> f64 { 1.616_255e-35 }

/// Planck mass [kg]
#[inline]
pub fn m_P() -> f64 { 2.176_434e-8 }

/// Planck time [s]
#[inline]
pub fn t_P() -> f64 { 5.391_247e-44 }

/// Convert electronvolts to joules
#[inline]
pub fn ev_to_j(ev: f64) -> f64 { ev * e() }

/// Convert joules to electronvolts
#[inline]
pub fn j_to_ev(j: f64) -> f64 { j / e() }

/// Convert femtometers to meters
#[inline]
pub fn fm_to_m(fm: f64) -> f64 { fm * 1e-15 }

/// Convert angstroms to meters
#[inline]
pub fn angstrom_to_m(a: f64) -> f64 { a * 1e-10 }
