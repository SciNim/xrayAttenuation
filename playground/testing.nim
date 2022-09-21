import std / [strformat]
import xrayAttenuation, ggplotnim


let h = Arsenic.init()
plotAttenuation(h)

Element[77].init().plotAttenuation()
echo readMolarMasses()

let ar = Argon.init()
let ρ_Ar = density(1050.mbar.to(Pascal), 293.K, ar.molarMass)
ar.plotAttenuation()
ar.plotTransmission(ρ_Ar, 3.cm.to(m))
echo absorptionLength(2.5.keV, numberDensity(ρ_Ar, ar.molarMass),
                      ar.f2eval(2.5.keV))

let au = Gold.init()
let ρ_Au = 19.32.g•cm⁻³
echo "num density Gold: ", numberDensity(ρ_Au, au.molarMass), " vs ", numberDensity(ρ_Ar, ar.molarMass)
#if true: quit()
let df = au.plotReflectivity(ρ_Au, 0.5.Degree) #, 0.5.Degree)

let Si₃N₄ = compound((Si, 3), (N, 4))
Si₃N₄.plotTransmission(3.44.g•cm⁻³, 300.nm.to(Meter))


block TestCompareWithHenke:
  let df = toDf({"Energy[keV]" : linspace(0.03, 10.0, 1000)})
    .mutate(f{float: "μ" ~ Si₃N₄.attenuationCoefficient(idx("Energy[keV]").keV).float},
            f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, 3.44.g•cm⁻³, 300.nm.to(Meter)).float},
            f{float: "Abs" ~ 1.0 - `Trans`})
  let dfR = readCsv("/home/basti/org/resources/Si3N4_density_3.44_thickness_0.3microns.txt", sep = ' ')
    .mutate(f{"Energy[keV]" ~ idx("PhotonEnergy(eV)") / 1000.0})
  echo dfR
  ggplot(df, aes("Energy[keV]")) +
    geom_line(aes = aes(y = "Trans")) +
    geom_line(data = dfR, aes = aes("Energy[keV]", "Transmission"), lineType = ltDashed, color = "red") +
    ggtitle(&"Transmission for Si₃N₄ thickness = {300.nm}, at ρ = {3.44.g•cm⁻³}") +
    ggsave("/tmp/transmission_si3n4_300nm.pdf", width = 800, height = 480)

block TestCompareAuReflectivityHenke:
  let dfR = readCsv(
    "../resources/reflectivity_gold_0_5deg.txt", sep = ' ', skipLines = 2,
    colNames = @["Energy[eV]", "Reflectivity"])
    .mutate(f{"Energy[keV]" ~ idx("Energy[eV]") / 1000.0})
  echo dfR
  ggplot(df, aes("Energy[keV]")) +
    geom_line(aes = aes(y = "Rs")) +
    xlim(0.0, 2.0) + ylim(0.85, 1.0) +
    geom_line(data = dfR, aes = aes("Energy[keV]", "Reflectivity"), lineType = ltDashed, color = "red") +
    ggtitle(&"Reflectivity Gold at θ = {0.5.°}, at ρ = {19.32.g•cm⁻³}") +
    ggsave("/tmp/reflectivity_gold_0_5.deg.pdf", width = 800, height = 480)
