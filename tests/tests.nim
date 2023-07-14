import std / [unittest, math]
import ../src/xrayAttenuation

import ggplotnim except almostEqual
import unchained

## A few simple tests

proc `=~=`[T: not seq](x, y: T): bool = almostEqual(x.float, y.float, epsilon = 3)
proc `=~=`[T](x, y: seq[T]): bool =
  check x.len == y.len
  result = true
  for i in 0 ..< x.len:
    result = result and x[i] =~= y[i]

suite "Basic physics utility procs":
  test "density":
    let Ar = Argon.init()
    #check density(1050.mbar.to(Pa), 293.15.K, Ar.molarMass) ==

  test "number density":
    let ρ = 1.g•cm⁻³
    #check numberDensity(ρ, M) ==

  test "wavelength":
    discard
    #check wavelength(1.keV) ==
    #check wavelength(5.keV) ==
    #check wavelength(10.keV) ==

  test "atomic absorption cross section":
    discard
    #check atomicAbsorptionCrossSection(1.keV,

  test "attenuation coefficient":
    #check attenuationCoefficient(1.keV, f2, M) =
    discard

  test "transmission":
    discard
    #check transmission(μ, ρ, length) ==

  test "absorption length":
    discard
    #check absorptionLength(energy, n_a, f2) ==

  test "delta":
    discard
    #check delta(energy, n_a, f1) ==

  test "beta":
    discard
    #check beta(energy, n_a, f1) ==

  test "refractive index":
    discard
    #check refractiveIndex*(energy, n_a, f0) ==

  test "wave number":
    discard
    #check waveNumber(energy, θ) ==

  test "reflectivity":
    discard
    # check reflectivity(θ, energy, n, σ) ==

  test "Depth graded multilayers":
    ## The following are the depth graded values for the LLNL telescope
    ## following a NuSTAR designed for the CAST experiment.
    check depthGradedLayers(d_min = 11.5.nm, d_max = 22.5.nm, N = 2, c = 1.0) =~= @[22.5.nm, 11.5.nm]
    check depthGradedLayers(d_min = 7.0.nm, d_max = 19.0.nm, N = 3, c = 1.0)  =~= @[19.nm, 10.2308.nm, 7.nm]
    check depthGradedLayers(d_min = 5.5.nm, d_max = 16.0.nm, N = 4, c = 1.0)  =~= @[16.nm, 9.77778.nm, 7.04.nm, 5.5.nm]
    check depthGradedLayers(d_min = 5.0.nm, d_max = 14.0.nm, N = 5, c = 1.0)  =~= @[14.nm, 9.65517.nm, 7.36842.nm, 5.95745.nm, 5.nm]

let ar = Argon.init()
var dfBoth = newDataFrame()
block PureArgon:
  let ρ_Ar = density(1050.mbar.to(Pascal), 293.K, ar.molarMass)
  #ar.plotAttenuation()
  #ar.plotTransmission(ρ_Ar, 3.cm.to(m))

  # compute absorption length for every energy and plot
  let df = toDf({ "Energy[keV]" : linspace(0.03, 10.0, 1000) })
    .mutate(f{float: "l_abs" ~ absorptionLength(idx("Energy[keV]").keV, numberDensity(ρ_Ar, ar.molarMass),
                                                ar.f2eval(idx("Energy[keV]").keV)).float },
            f{"Type" <- "Argon"})
  dfBoth.add df
  ggplot(df, aes("Energy[keV]", "l_abs")) +
    geom_line() +
    ylab("l_abs [m]") +
    scale_y_log10() +
    ggtitle("Absorption length in Argon at 1050 mbar") +
    ggsave("/tmp/absorption_length_argon_1050mbar.pdf")


block ArAsGasMixture:
  let arC = compound((Ar, 1)) # need Argon gas as a Compound
  let gm = initGasMixture(293.K, 1050.mbar, [(arC, 1.0)])
  let df = toDf({ "Energy[keV]" : linspace(0.03, 10.0, 1000) })
      .mutate(f{float: "l_abs" ~ absorptionLength(gm, idx("Energy[keV]").keV).float})
  ggplot(df, aes("Energy[keV]", "l_abs")) +
    geom_line() +
    ylab("l_abs [m]") +
    scale_y_log10() +
    ggtitle("Absorption length in Argon at 1050 mbar") +
    ggsave("/tmp/absorption_length_argon_as_mixt_1050mbar.pdf")



block ArgonIsobutane:
  let arC = compound((Ar, 1)) # need Argon gas as a Compound
  let isobutane = compound((C, 4), (H, 10))
  let gm = initGasMixture(293.K, 1050.mbar, [(arC, 0.977), (isobutane, 0.023)])

  # compute absorption length for every energy and plot
  let df = toDf({ "Energy[keV]" : linspace(0.03, 10.0, 1000) })
    .mutate(f{float: "l_abs" ~ absorptionLength(gm, idx("Energy[keV]").keV).float},
            f{"Type" <- "Argon/Isobutane"})
  dfBoth.add df
  echo df
  ggplot(df, aes("Energy[keV]", "l_abs")) +
    geom_line() +
    ylab("l_abs [m]") +
    scale_y_log10() +
    ggtitle("Absorption length in Argon/Isobutane at 1050 mbar") +
    ggsave("/tmp/absorption_length_argon_isobutane_1050mbar.pdf")



ggplot(dfBoth, aes("Energy[keV]", "l_abs", color = "Type")) +
  geom_line() +
  ylab("l_abs [m]") +
  scale_y_log10() +
  ggtitle("Absorption length in Argon (/Isobutane) at 1050 mbar") +
  ggsave("/tmp/absorption_length_argon_and_ar_isobutane_1050mbar.pdf")

echo getFluorescenceLines(ar)
echo getFluorescenceLines(Silver.init())
