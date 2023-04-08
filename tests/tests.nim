# import unittest

import xrayAttenuation, ggplotnim, unchained


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
