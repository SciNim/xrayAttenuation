import unchained, math, complex

## This file implements helpers as well as all the base calculations
## required to compute transmissions / absorptions / reflections. The
## main file implements a convenience layer on top of these that take
## `Element/Compound` types to automate the majority of the argument
## handling.
##
## Note: the naming of the procedures in here is verbose for clarity.
## More common abbreviations may be added as templates in the future.

let R* = N_A * k_B # 8.314.J•K⁻¹•mol⁻¹ the universal Gas constant

## A bunch of units we use in this library. They have to be defined here,
## as they are used as arguments to procedures, which requires them to be
## predefined.
#defUnit(m², toExport = true)
defUnit(cm⁻¹, toExport = true)
defUnit(cm², toExport = true)
defUnit(m⁻², toExport = true)
defUnit(cm⁻³, toExport = true)
defUnit(cm²•g⁻¹, toExport = true)
defUnit(g•cm⁻³, toExport = true)
defUnit(g•mol⁻¹, toExport = true)
defUnit(m⁻³, toExport = true)
defUnit(m⁻³•kg, toExport = true)
defUnit(m⁻³•mol, toExport = true)
defUnit(m⁻¹, toExport = true)
defUnit(Joule•Meter, toExport = true)

proc density*[P: Pressure](p: P, T: Kelvin, M: g•mol⁻¹): g•cm⁻³ =
  ## Compute the density `ρ` using the ideal gas law based on the given
  ## pressure `P` (use `.to(Pascal)` in the argument if you have a different
  ## unit), temperature `T` and molar mass `M`.
  result = (p * M / (R * T)).to(g•cm⁻³)

proc numberDensity*(ρ: g•cm⁻³, M: g•mol⁻¹): cm⁻³ =
  ## Computes the number density of a medium with density `ρ` and molar mass
  ## `M`. That is the number of particles in a unit volume (of `cm⁻³` in this case).
  result = N_A / M * ρ

proc wavelength*(energy: keV): Meter =
  ## Compute the wavelength of the X-ray with the given `energy`.
  result = (hp * c / energy).to(Meter)

proc atomicAbsorptionCrossSection*(energy: keV, f2: UnitLess): cm² =
  ## Computes the atomic absoprtion cross section `σ_a` based on the scattering factor `f2`
  ## via
  ##  `σ_A = 2 r_e λ f₂`
  ## for the given `energy`.
  let λ = wavelength(energy)
  result = (2 * r_e * λ * f2).to(cm²)

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
  let σ_a = atomicAbsorptionCrossSection(energy, f2)
  result = (σ_a * N_A / M)

proc transmission*(μ: cm²•g⁻¹, ρ: g•cm⁻³, length: Meter): UnitLess =
  ## Computes the transmission given a mass absorption coefficient `μ`, a density
  ## of the medium `ρ` and some `length` according to the Beer-Lambert law.
  ##
  ## Note: again, take care to convert the length dimensions correctly!
  result = exp(-μ * ρ * length)

proc absorptionLength*(energy: keV, n_a: cm⁻³, f2: UnitLess): Meter =
  ## Computes the absorption length given a number density `n_a` and
  ## scattering factor `f2` at an energy `keV`.
  ##
  ## Equivalent to `1.0 / (μ * ρ)`.
  let σ_A = atomicAbsorptionCrossSection(energy, f2)
  result = (1.0 / (n_a * σ_A)).to(Meter)

proc delta*(energy: keV, n_a: cm⁻³, f1: float): float =
  ## Computes `delta` at the given `energy` given scattering factor `f1` and
  ## number density `n_a`.
  ##
  ## `δ` is the real portion of the (1 - refractive index) formulation
  ## and defined by:
  ##   `δ = n_a r_e λ² / 2π f_1(ω)`
  let λ = wavelength(energy)
  result = n_a * r_e * λ*λ / (2 * π) * f1

proc beta*(energy: keV, n_a: cm⁻³, f2: float): float =
  ## Computes `beta` at the given `energy` given scattering factor `f2` and
  ## number density `n_a`.
  ##
  ## `β` is the imaginary portion of the (1 - refractive index) formulation
  ## and defined by:
  ##
  ## `β = n_a r_e λ² / 2π f_2(ω)`
  let λ = wavelength(energy)
  result = n_a * r_e * λ*λ / (2 * π) * f2

proc refractiveIndex*(energy: keV, n_a: cm⁻³, f0: Complex[float]): Complex[float] =
  ## Computes the complex refractive index for the given `energy` based on the
  ## the number density of the material `n_a` and combined scattering factor
  ##
  ##  `f0(ω) = f1(ω) + if2(ω)`
  ##
  ##  `n(ω) = 1 - r_e λ² / 2π Σ_i n_ai fi(0)`
  ## where the sum runs over all constituents of the material with `fi(0)` the
  ## complex scattering factor (forward scattering) and `n_ai` the number density
  ## of the constituent element `i`.
  ##
  ## NOTE: currently only for single elements
  let λ = wavelength(energy)
  result = 1.0 - (r_e * λ * λ / (2 * π) * n_a).float * f0

proc waveNumber*(energy: keV, θ: Degree): m⁻¹ =
  ## Computes the wave number `k` for an incoming wave with `energy`
  ## and incident angle `θ` (measured from surface).
  result = 2*π * sin(θ.to(Radian)) / wavelength(energy)

proc reflectivity*(θ: Degree, energy: keV, n: Complex[float], σ: Meter): float =
  ## Computes the reflectivty `R = |r²|` for the X-ray at the given `energy`
  ##
  ## `σ`: surface roughness
  ##
  ## TODO: understand why other equation doesn't work, check that surface roughness
  ## correction works as expected (and implement it using this approach)
  ## `R = |r_0 exp(-2 k_i k_z σ²)|²`
  let k = waveNumber(energy, θ).float
  let kθ = k * cos(θ.to(Radian))

  let km = sqrt( k*k - kθ*kθ )
  let kp = sqrt( (k*k).float * n*n - kθ*kθ )
  result = abs2( (km - kp) / (km + kp) )

proc rij_s*(ni, nj: float, θ: Degree): float =
  let niS = ni * sin(θ.to(Radian))
  let njS = nj * sin(θ.to(Radian))
  result = (niS - njS) / (niS + njS)

#proc refR*(

proc scatteringPotential(ρ: cm⁻³): m⁻² =
  result = (4*π * r_e * ρ).to(Meter⁻²)

proc reflectivity*(θ: Degree, ρ: cm⁻³, energy: keV, n: Complex[float], σ: Meter): float =
  ## `σ`: surface roughness
  let λ = wavelength(energy)

  let sθ = sin(θ.to(Radian))
  let sθp = sqrt( sθ*sθ - (1.0 - n)^2 ) / (n^2)

  let k = waveNumber(energy, θ).float
  let kθ = k * cos(θ.to(Radian))
  #let ks = n * k * sθp
  #let ks = sqrt( k*k - scatteringPotential(ρ).float )#k*k * sθ*sθ - k*k * (1.0 - n)*(1.0 - n) )
  #result = abs( (k - ks) / (k + ks) )^2
  #result = abs2((sθ - n * sθp) / (sθ + n * sθp))


  ## XXX: THIS FINALLY WORKS. Well, almost. The numbers are a bit too small.
  # figure out what's wrong with the other ones. Take pen and paper.
  let km = sqrt( k*k - kθ*kθ )
  let kp = sqrt( (k*k).float * n*n - kθ*kθ )
  result = abs2( (km - kp) / (km + kp) )  # dumb
  when false:
    let cθ = cos(θ.to(Radian))
    let k_iz = (2*π / λ * cθ).float
    let k_tz = (2*π / λ).float * sqrt( n*n - cθ*cθ ) #sqrt( n*n - cθ*cθ ) #sqrt( n*n - cθ*cθ ) #complex(- cθ*cθ, 0.0) )
    echo "kiz = ", kiz, " ktz = ", ktz, " from ", abs2(n), " c^2(θ) = ", cθ * cθ
    let r0 = abs2( (k_iz - k_tz) / (k_iz + k_tz) )
    result = r0 # * r0 #abs2(r0 * exp(- 2 * k_iz * k_tz * (σ^2).float ))
