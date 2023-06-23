import std / [strformat]
import xrayAttenuation, ggplotnim


let h = Arsenic.init()
plotAttenuation(h)

Element[77].init().plotAttenuation()
echo readMolarMasses()

let ar = Argon.init()
let ρ_Ar = density(1050.mbar.to(Pascal), 293.K, ar.molarMass)
#let ρ_Ar = density((1400 * 0.973).mbar.to(Pascal), 300.K, ar.molarMass)
ar.plotAttenuation(logLog = true)
## XXX:: the following with typo `enrgyMin` causes a compiler segfault!
## ar.plotTransmission(ρ_Ar, 3.cm.to(m), enrgyMin = 5.5, energyMax = 6.5)
ar.plotTransmission(ρ_Ar, 3.cm.to(m), energyMin = 5.5, energyMax = 6.5)
#block ArgonTrans:
#  let μ_Ar = ar.attenuationCoefficient(5.9.keV)
#  let t = transmission(μ_Ar, ρ_Ar, 3.cm.to(Meter))
#  echo "Trasmission at 5.9 keV = ", t
#
#echo absorptionLength(2.5.keV, numberDensity(ρ_Ar, ar.molarMass),
#                      ar.f2eval(2.5.keV))
#echo ar.absorptionLength(ρ_Ar, 2.5.keV)
#
let au = Gold.init()
let ρ_Au = 19.32.g•cm⁻³
echo "num density Gold: ", numberDensity(ρ_Au, au.molarMass), " vs ", numberDensity(ρ_Ar, ar.molarMass)
##if true: quit()
let df = au.plotReflectivity(ρ_Au, 0.5.Degree) #, 0.5.Degree)
#
#let Si₃N₄ = compound((Si, 3), (N, 4))
#Si₃N₄.plotTransmission(3.44.g•cm⁻³, 300.nm.to(Meter))
#
let mylar = compound((C, 10), (H, 8), (O, 4))
mylar.plotTransmission(1.4.g•cm⁻³, 4.μm.to(Meter), energyMax = 10.0)
#
#echo Si₃N₄.absorptionLength(3.44.g•cm⁻³, 2.5.keV)
#
#block TestCompareWithHenke:
#  let df = toDf({"Energy[keV]" : linspace(0.03, 10.0, 1000)})
#    .mutate(f{float: "μ" ~ Si₃N₄.attenuationCoefficient(idx("Energy[keV]").keV).float},
#            f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, 3.44.g•cm⁻³, 300.nm.to(Meter)).float},
#            f{float: "Abs" ~ 1.0 - `Trans`})
#  let dfR = readCsv("/home/basti/org/resources/Si3N4_density_3.44_thickness_0.3microns.txt", sep = ' ')
#    .mutate(f{"Energy[keV]" ~ idx("PhotonEnergy(eV)") / 1000.0})
#  echo dfR
#  ggplot(df, aes("Energy[keV]")) +
#    geom_line(aes = aes(y = "Trans")) +
#    geom_line(data = dfR, aes = aes("Energy[keV]", "Transmission"), lineType = ltDashed, color = "red") +
#    ggtitle(&"Transmission for Si₃N₄ thickness = {300.nm}, at ρ = {3.44.g•cm⁻³}") +
#    ggsave("/tmp/transmission_si3n4_300nm.pdf", width = 800, height = 480)
#
#block TestCompareAuReflectivityHenke:
#  let dfR = readCsv(
#    "../resources/reflectivity_gold_0_5deg.txt", sep = ' ', skipLines = 2,
#    colNames = @["Energy[eV]", "Reflectivity"])
#    .mutate(f{"Energy[keV]" ~ idx("Energy[eV]") / 1000.0})
#  echo dfR
#  ggplot(df, aes("Energy[keV]")) +
#    geom_line(aes = aes(y = "Rs")) +
#    xlim(0.0, 2.0) + ylim(0.85, 1.0) +
#    geom_line(data = dfR, aes = aes("Energy[keV]", "Reflectivity"), lineType = ltDashed, color = "red") +
#    ggtitle(&"Reflectivity Gold at θ = {0.5.°}, at ρ = {19.32.g•cm⁻³}") +
#    ggsave("/tmp/reflectivity_gold_0_5.deg.pdf", width = 800, height = 480)


import sequtils, ggplotnim
block LLNL_Reflectivity:
#when false:
  #[
Reference:
let m1 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 11.5,
                        D_max = 22.5,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=2,
                        SigmaValues=[1.0])
]#

  # initalize the elements of the substrate
  let C = Carbon.init()
  let Pt = Platinum.init()
  let glass = compound((Si, 1), (O, 2)) # pure quartz glass
  const ρC = 2.26.g•cm⁻³ # 2.267 should be correct. used number same as DarpanX however
  const ρPt = 21.41.g•cm⁻³ # 21.45 should be correct. used number same as DarpanX however
  const ρGlass = 2.65.g•cm⁻³
  const θ = 0.5.°

  if false:
    echo 0.03.keV.wavelength().m.to(nm)
    echo "to"
    echo 10.keV.wavelength().m.to(nm)
    quit()

  block RefractiveIndexDebug:
    let λ = 1.976788861549859121e2.nm
    echo λ.energy()
    let nC = refractiveIndex(C, λ.energy(), ρC)
    echo "Ours = ", nC
    let target = complex(1.705354170365389610, 0.000000000000000000)
    echo "Target = ", target

  block RefrIndexGlass:
    echo glass.molarWeight(), " to target: #A(g/mol)=60.0843"
    #[
#Formula=SiO2
#A(g/mol)=60.0843
#Rho(g/cm3)=2.65
#Energy range=0.006958-432.9451 keV
#
# Wavelength (A) 	 n  	 k
1.781894185059028132e+03 2.009538984318353094e+00 5.633583284510720324e-01
1.755031458651102866e+03 2.740088255913744497e+00 5.342169546214327047e-01
1.748004305663511332e+03 3.617751372225351147e+00 5.267290442116252525e-01
]#
    # let λ = 1.781894185059028132e2.nm

    proc checkLambda(λ: nm, target: Complex[float]) =
      echo "At energy: ", λ.energy()
      let nGlass = refractiveIndex(glass, λ.energy(), ρGlass)
      echo "⇒⇒⇒⇒Ours = ", nGlass
      echo "⇒⇒⇒⇒Target = ", target
    let λ = 1.516483801844758261.nm
    let target = complex(9.991328857359060844e-1, 1.790930224679519115e-4) #complex(2.009538984318353094, 5.633583284510720324e-01)
    checkLambda(λ, target)
    let λ2 = 4.263076242678597509e1.nm
    let t2 = complex(7.746077971982764376e-01, 1.683285257772849097e-01)
    let λ3 = 3.987909862920872683e+01.nm
    let t3 = complex(7.946499269987641334e-01, 1.483447463676793865e-01)
    checkLambda(λ2, t2)
    checkLambda(λ3, t3)
    let λ4 = 2.856648415496417215e+01.nm
    let t4 = complex(8.737283973357047318e-01, 7.454228279449584549e-02)
    checkLambda(λ4, t4)


  #if true: quit()

  let Γ = 0.45
  let dMin = 11.5.nm
  let dMax = 22.5.nm
  let N = 2
  var ds = @[dMax * Γ, dMax * (1.0 - Γ), dMin * Γ, dMin * (1.0 - Γ)]
  let nVacuum = complex(1.0, 0.0) # vacuum
  let energies = linspace(0.03, 10.0, 1000).mapIt(it.keV)
  var RMs = newSeq[float]()
  for i, E in energies:
    # 1. compute the refractive indices for each layer
    let nC = refractiveIndex(C, E, ρC)
    let nPt = refractiveIndex(Pt, E, ρPt)
    let nGlass = refractiveIndex(glass, E, ρGlass)
    echo "-------------------- \n\n\n\n nC ", nC, " nPt ", nPt, " nGlass ", nGlass
    let ns = @[nVacuum, nC, nPt, nC, nPt, nGlass]
    RMs.add multilayerReflectivity(θ, E, ns, ds, parallel = false)

  ggplot(toDf({"E" : energies.mapIt(it.float), "R" : RMs}), aes("E", "R")) +
    geom_line() +
    ggtitle("Reflectivity of LLNL recipe 1") +
    ggsave("/tmp/llnl_reflectivity_1.pdf")



#block NustarReflectivity:
#  # initalize the elements of the substrate
#  # m=drp.Multilayer(MultilayerType="DepthGraded",SubstrateMaterial="SiO2",LayerMaterial=["Pt","C"],Repetition=150,D_max=128.1,D_min=31.7,C=0.245,GammaTop=0.7,Gamma=0.45,SigmaValues=[4.5])
#  let C = Carbon.init()
#  let Pt = Platinum.init()
#  let glass = compound((Si, 1), (O, 2)) # pure quartz glass
#  const ρC = 2.26.g•cm⁻³ # 2.267 should be correct. used number same as DarpanX however
#  const ρPt = 21.41.g•cm⁻³ # 21.45 should be correct. used number same as DarpanX however
#  const ρGlass = 2.65.g•cm⁻³
#  const θ = 0.5.°
#
#  let Γ = 0.45
#  let dMin = 11.5.nm
#  let dMax = 22.5.nm
#  let N = 2
#  var ds = @[10000.nm, dMax * Γ, dMax * (1.0 - Γ), dMin * Γ, dMin * (1.0 - Γ)]
#  let nVacuum = complex(1.0, 0.0) # vacuum
#  let energies = linspace(0.03, 10.0, 1000).mapIt(it.keV)
#  var RMs = newSeq[float]()
#  for i, E in energies:
#    # 1. compute the refractive indices for each layer
#    let nC = refractiveIndex(C, E, ρC)
#    let nPt = refractiveIndex(Pt, E, ρPt)
#    let nGlass = refractiveIndex(glass, E, ρGlass)
#    echo "-------------------- \n\n\n\n nC ", nC, " nPt ", nPt, " nGlass ", nGlass
#    let ns = @[nVacuum, nC, nPt, nC, nPt, nGlass]
#    RMs.add multilayerReflectivity(θ, E, ns, ds)
#
#  ggplot(toDf({"E" : energies.mapIt(it.float), "R" : RMs}), aes("E", "R")) +
#    geom_line() +
#    ggtitle("Reflectivity of LLNL recipe 1") +
#    ggsave("/tmp/llnl_reflectivity_1.pdf")
