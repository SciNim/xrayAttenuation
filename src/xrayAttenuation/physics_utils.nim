import unchained, math

let R* = N_A * k_B #8.314.J•K⁻¹•mol⁻¹

defUnit(cm⁻³, toExport = true)
defUnit(cm²•g⁻¹, toExport = true)
defUnit(g•cm⁻³, toExport = true)
defUnit(g•mol⁻¹, toExport = true)
proc density*(p: Pascal, T: Kelvin, M: g•mol⁻¹): g•cm⁻³ =
  result = (p * M / (R * T)).to(g•cm⁻³)

proc wavelength*(energy: keV): Meter =
  result = hp * c / (energy.to(Joule))

defUnit(m⁻³, toExport = true)
defUnit(m⁻³•kg, toExport = true)
defUnit(m⁻³•mol, toExport = true)
defUnit(m⁻¹, toExport = true)
proc attenuationCoefficient*(energy: keV, f2: UnitLess, M: g•mol⁻¹): cm²•g⁻¹ =
  ## Computes the mass attenuation coefficient μ_m based on the scattering factor `f2`
  ## via
  ##  `σ_A = 2 r_e λ f₂`
  ## where `r_e` is the classical electron radius, `λ` the wavelength of the photon
  ## considered and `f₂` the scattering factor of the material (a unitless quantity).
  ## `σ_A` then is the atomic photoabsorption cross section.
  ##
  ## Further,
  ##  `μ_m = N_A / M * σ_A`
  ## then defines the mass absorption coefficient in units of `cm²•g⁻¹` (as the molar
  ## mass `M` has units `g•mol⁻¹`).
  ##
  ## Note: If implemented from hand, be careful with the units of `r_e` and `λ` (don't
  ## forget to convert from `m` to `cm`!
  ##
  ## Note 2: The Henke notes here:
  ## https://xdb.lbl.gov/Section1/Sec_1-6.pdf
  ## use `A` to refer to `M` (`M = A * M_u`, where `M_u` is the molar mass constant), which
  ## is extremely confusing as in a related document:
  ## https://xdb.lbl.gov/Section5/Sec_5-5.pdf (eq. 3.26)
  ## they refer to the real `A` (i.e. unit less standard atomic weight) using `A`!
  result = (2 * r_e * wavelength(energy) * f_2 * N_A / M).to(cm²•g⁻¹)

proc transmission*(μ: cm²•g⁻¹, ρ: g•cm⁻³, length: Meter): UnitLess =
  ## Computes the transmission given a mass absorption coefficient `μ`, a density
  ## of the medium `ρ` and some `length` according to the Beer-Lambert law.
  ##
  ## Note: again, take care to convert the length dimensions correctly!
  result = exp(-μ * ρ * length)

proc numberDensity*(ρ: g•cm⁻³, M: g•mol⁻¹): cm⁻³ =
  ## Computes the number density of a medium with density `ρ` and molar mass
  ## `M`. That is the number of particles in a unit volume (of `cm⁻³` in this case).
  result = N_A / M * ρ

proc absorptionLength*(energy: keV, n_a: cm⁻³, f2: UnitLess): Meter =
  ## Computes the absorption length given a number density `n_a` and
  ## scattering form factor `f2` at an energy `keV`.
  result = 1.0 / (2 * n_a * r_e * wavelength(energy) * f2)
