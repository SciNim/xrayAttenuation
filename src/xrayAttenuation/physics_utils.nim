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

proc energy*[L: Length](λ: L): keV =
  ## Compute the wavelength of the X-ray with the given `energy`.
  result = (hp * c / λ).to(keV)

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
  result = 1.0 - (n_a * r_e * λ*λ / (2 * π) ) * f0

proc waveNumber*(energy: keV, θ: Degree): m⁻¹ =
  ## Computes the wave number `k` for an incoming wave with `energy`
  ## and incident angle `θ` (measured from surface).
  result = 2*π * sin(θ.to(Radian)) / wavelength(energy)

proc refractedAngle*(θi: Degree, n_i, n_j: Complex[float]): Complex[float] =
  ## Given the incidence angle `θi` returns the refracted angle according
  ## to Snell's law for the materials with refractive indices `n_i` and `n_j`
  ##
  ## `n_i sin(θ_i) = n_j sin(θ_j)`
  ## ⇒ `θ_j = arcsin(n_i sin(θ_i) / n_j)`
  ##
  ## `IMPORTANT`: The incident angle in this formulation is measured from the
  ## normal of the interface! For grazing angles we typically measure from the
  ## interface, so convert `90.° - θi` before calling!
  ##
  ## `IMPORTANT 2`: This procedure is unsafe in the contexts of total reflection,
  ## i.e. for grazing angles if the `n_i > n_j` (vacuum to metal in X-ray regime
  ## or glass to vacuum in visible light for example).
  result = arcsin(n_i * sin(θ_i.to(Radian)) / n_j)

proc refractedAngleSin*(θi: Degree, n_i, n_j: Complex[float]): Complex[float] =
  ## Given the incidence angle `θi` returns the sine of the refracted angle according
  ## to Snell's law for the materials with refractive indices `n_i` and `n_j`
  ##
  ## `n_i sin(θ_i) = n_j sin(θ_j)`
  ## ⇒ `sin(θ_j) = n_i sin(θ_i) / n_j`
  ##
  ## `IMPORTANT`: The incident angle in this formulation is measured from the
  ## normal of the interface! For grazing angles we typically measure from the
  ## interface, so convert `90.° - θi` before calling!
  result = n_i * sin(θ_i.to(Radian)) / n_j

proc refractedAngleSin*(sinθi: Complex[float], n_i, n_j: Complex[float]): Complex[float] =
  ## Version of the above if `sin(θ_i)` is already given as complex number.
  result = n_i * sinθi / n_j

proc scatteringPotential(ρ: cm⁻³): m⁻² =
  result = (4*π * r_e * ρ).to(Meter⁻²)

proc reflectivity*(sinθn_i, n_i, sinθn_j, n_j: Complex[float], energy: keV, σ: Meter,
                   parallel: bool): Complex[float] =
  ## Computes the actual reflectvity given incidence and refracted angles as the
  ## `sin(θ)` and the refractive indices. The angles are given from the normal
  ## of the interface.
  ##
  ## If `parallel` is `true` computes the reflectivtiy for `p` polarization, else
  ## for `s` pol.
  ##
  ## `σ` is the surface roughness.
  ##
  ## TODO: IMPLEMENT `σ` correction!
  ##
  ## The calculation is done using the Fresnell equations, which relate the
  ## refractive indices and the incident / outgoing angles to the reflectivity
  ## and transmission.
  ## To use these equations for grazing angles for X-rays (in which the real part
  ## of the refracted angle is often < 1), it is important to keep the angles
  ## as the `sin(θ)` to avoid having to compute `arcsin( <expr> )`, which may
  ## be undefined.
  ##
  ## For a directly applicable treatment of the math required for this, see
  ##
  ## `David L. Windt, 1998 - IMD - Software for modeling the optical properties of multilayer films`
  ##
  ## Fresnell equations:
  ##
  ## s-polarization:
  ##
  ## `         n_i · cos(θ_i) - n_j · cos(θ_j)`
  ## `r^s_ij = ------------------------------ `
  ## `         n_i · cos(θ_i) + n_j · cos(θ_j)`
  ##
  ## p-polarization:
  ##
  ## `         n_i · cos(θ_j) - n_j · cos(θ_i)`
  ## `r^p_ij = ------------------------------ `
  ## `         n_i · cos(θ_j) + n_j · cos(θ_i)`
  ##
  ## where `θ_i`, `θ_j` are the incident, refracted angles and `n_i`, `n_j` the
  ## refractive indices on the incident / outgoing side.
  ##
  ## Given that we have to use the `sin(θ)` to express angles, we use the identity
  ##
  ## `sin²(x) + cos²(x) = 1`
  ##
  ## to express the `cos(θ)` as `√(1 - sin²(θ))` for the equation.
  let cθi = sqrt(1.0 - sinθn_i*sinθn_i)
  let cθj = sqrt(1.0 - sinθn_j*sinθn_j)
  if parallel:
    result = (n_i * cθj - n_j * cθi) / (n_i * cθj + n_j * cθi)
  else:
    result = (n_i * cθi - n_j * cθj) / (n_i * cθi + n_j * cθj)

proc reflectivity*(θ_i: Degree, energy: keV, n_j: Complex[float], σ: Meter,
                   parallel: bool): float =
  ## Computes the reflectivity of the grazing angle `θ_i` (from the interface) given
  ## the medium `n_j`. The medium on the incident side is assumed to be vacuum.
  ##
  ## To compute the reflectivity more generally with other media, use the
  ## version above which takes the sine of angles from the normal of the interface.
  ##
  ## If `parallel` is `true` computes the reflectivtiy for `p` polarization, else
  ## for `s` pol.
  ##
  ## `σ` is the surface roughness.
  ##
  ## (TODO: perform surface roughness correction via
  ##
  ## `r'_ij = r_ij * exp(- s_i² σ² / 2)`
  ##
  ## with `s_i`
  ##
  ## `s_i = 4π cos(θ_i) / λ`
  let θn_i = 90.° - θ_i # incidence angle to normal vector
  let sinθn_i = complex(sin(θn_i.to(Radian)), 0.0)
  let n_i = complex(1.0, 0.0) ## Vacuum refractive index
  let sinθn_j = refractedAngleSin(sinθn_i, n_i, n_j)
  let r = reflectivity(sinθn_i, n_i, sinθn_j, n_j, energy, σ, parallel)
  result = abs2(r)

proc reflectivity*(θ_i: Degree, energy: keV, n_j: Complex[float], σ: Meter): float =
  ## Overload of the above, which assumes the incoming light is unpolarized, i.e.
  ## returns the average of the s- and p-polarizations:
  ##
  ## `R = 1/2 (|r_s|² + |r_p|²)`
  let r_s = reflectivity(θ_i, energy, n_j, σ, parallel = false)
  let r_p = reflectivity(θ_i, energy, n_j, σ, parallel = true)
  result = 0.5 * (r_s + r_p)

proc reflectivityAlt*(θ: Degree, energy: keV, n: Complex[float], σ: Meter): float =
  ## Alternative implementation via the wavenumber of the reflectivity.
  ##
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

proc multilayerReflectivity*(θ_i: Degree, energy: keV, ns: seq[Complex[float]], ds: seq[nm],
                             parallel: bool): float =
  ## Incoming angle from interface `θ_i`.
  ##
  ## Multilayer consisting of the refractive indices given by `ns` where
  ## `ns[0]` is the ambient medium (likely air or vacuum), `ns[^1]` is the
  ## substrate medium and the rest are the actual interfaces that should
  ## produce reflection.
  ## `ds` are the thicknesses of the different layers. Note that
  ## `ds.len == ns.len - 2`.
  ##
  ## Computes it recursively by
  ## ```
  ##        r_ij + r_j exp(2 i β_i)
  ## r_i = -------------------------
  ##       1 + r_ij r_j exp(2 i β_i)
  ## ```
  ## with
  ##
  ## `β_i = 2*π * d_i * cos(θ_i) / λ` (`θ_i` from normal axis!)
  ##
  ## and
  ##
  ## r_ij computed via `reflectivity`
  ## (TODO: perform surface roughness correction in `reflectivity`)
  ##
  ## `r_i` is computed for s- and p polarizations separately!

  if ns.len < 3:
    raise newException(ValueError, "Need at least 3 input refraction indices. One for the ambient, " &
      "one for the substrate and one or more for the layers to reflect on.")
  if ds.len != ns.len - 2:
    raise newException(ValueError, "Need thicknesess for each layer, i.e. 2 less than number of " &
      "refractive indices. Got: " & $ds & " for number of refractive indices: " & $ns.len)

  proc beta(sinθ_i, n_i: Complex[float], d_i, λ: nm): Complex[float] =
    ## Given layer thickness `d_i` and wavelength `λ` and incidence angle
    ## `θ_i` (from normal and not from interface!), compute the parameter
    ## `β_i`.
    ##
    ## Note: merged factor `2` from `exp(2 i β_i)` into definition of `β_i`
    let cθi = sqrt(1.0 - sinθ_i*sinθ_i)
    let d_by_λ = d_i / λ
    result = 4*π * d_by_λ * n_i * cθi

  # 1. convert incident angle to normal of interface and take `sin`
  ## NOTE: angles are actually the sin(θ) and a complex number! This is to deal with cases
  ## correctly in which grazing angle so small as the produce total reflection where Snell's law
  ## is not valid for real numbers.
  var θn_i = complex(sin((90.° - θ_i).to(Radian)), 0.0) # incidence to normal angle, start from incidence to material
  let sθn_0 = θn_i

  let λ = wavelength(energy)
  ## NOTE: We can compute the refracted angle for any layer based on the 0-th layer, because
  ## `n_i · sin(θ_i) = n_j · sin(θ_j) = n_k · sin(θ_k)`
  # 2. compute incident angle to substrate, i.e.
  # = equal to refracted angle of layer _before_ substrate
  let sθ_last = refractedAngleSin(sθn_0, ns[0], ns[^2])
  # 3. refracted angle in substrate
  let sθ_substrate = refractedAngleSin(sθn_0, ns[0], ns[^1])
  # 4. compute `r_j` (= reflectivity) of the substrate
  var r_j = reflectivity(sθ_last, ns[^2], sθ_substrate, ns[^1], energy, 0.0.m, parallel = parallel)
  # 5. loop back from layer _above_ substrate to first layer. `ns` has `N + 2`
  # elements (ns[0] = vacuum / ambient, `ns[^1]` substrate), i.e. loop all medium
  # layers from bottom up.
  for i in countdown(ns.high - 1, 1):
    # 1. compute `β_i`. Need incidence angle of current layer `i` and refractive index
    let sθi = ns[0] * sθn_0 / ns[i]
    let β_i = beta(sθi, ns[i], ds[i-1], λ.to(NanoMeter))
    # 2. compute `r_ij`
    let angl = ns[0] * sθn_0 / ns[i-1]
    let r_ij = reflectivity(angl, ns[i-1], sθi, ns[i], energy, 0.0.m, parallel = parallel)
    # 3. assemble `r_j`
    r_j = (r_ij + r_j * exp(im(1.0) * β_i)) / (1.0 + r_ij * r_j * exp(im(1.0) * β_i))
  # 4. once we are at the end of the loop, the last `r_j` is our final reflectivity
  result = abs2(r_j)
