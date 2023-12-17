import std / [strformat, math]
import xrayAttenuation, ggplotnim

import sequtils, ggplotnim
#when false:
  #[
Reference:

let m1 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 115.0,
                        D_max = 225.0,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=2,
                        SigmaValues=[1.0])

let m2 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 70.0,
                        D_max = 190.0,
                        Gamma = 0.45,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=3,
                        SigmaValues=[1.0])

let m3 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 55.0,
                        D_max = 160.0,
                        Gamma = 0.4,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=4,
                        SigmaValues=[1.0])

let m4 = drp.Multilayer(MultilayerType="DepthGraded",
                        SubstrateMaterial="SiO2",
                        D_min = 50.0,
                        D_max = 140.0,
                        Gamma = 0.4,
                        C = 1.0,
                        LayerMaterial=["C", "Pt"],
                        Repetition=5,
                        SigmaValues=[1.0])

]#

proc manual() =
  const ρC = 2.26.g•cm⁻³ # 2.267 should be correct. used number same as DarpanX however
  const ρPt = 21.41.g•cm⁻³ # 21.45 should be correct. used number same as DarpanX however
  const ρGlass = 2.65.g•cm⁻³

  let C = Carbon.init(ρC)
  let Pt = Platinum.init(ρPt)
  let glass = compound((Si, 1), (O, 2), ρ = ρGlass) # pure quartz glass
  const θ = 0.5.°

  let Γ = 0.45
  let dMin = 11.5.nm
  let dMax = 22.5.nm
  let N = 2
  var ds = @[dMax * Γ, dMax * (1.0 - Γ), dMin * Γ, dMin * (1.0 - Γ)]
  let nVacuum = complex(1.0, 0.0) # vacuum
  let energies = linspace(0.03, 10.0, 1000).mapIt(it.keV)
  var RMs = newSeq[float]()
  for i, E in @[1.5.keV]: #energies:
    # 1. compute the refractive indices for each layer
    let nC = refractiveIndex(C, E, ρC)
    let nPt = refractiveIndex(Pt, E, ρPt)
    let nGlass = refractiveIndex(glass, E, ρGlass)
    echo "-------------------- \n\n\n\n nC ", nC, " nPt ", nPt, " nGlass ", nGlass
    let ns = @[nVacuum, nC, nPt, nC, nPt, nGlass]
    RMs.add multilayerReflectivity(θ, E, ns, ds, parallel = false)

  #ggplot(toDf({"E" : energies.mapIt(it.float), "R" : RMs}), aes("E", "R")) +
  #  geom_line() +
  #  ggtitle("Reflectivity of LLNL recipe 1") +
  #  ggsave("/tmp/llnl_reflectivity_1.pdf")


template getLayers(): untyped =
  const ρC = 2.26.g•cm⁻³ # 2.267 should be correct. used number same as DarpanX however
  const ρPt = 21.41.g•cm⁻³ # 21.45 should be correct. used number same as DarpanX however
  const ρGlass = 2.65.g•cm⁻³

  let C = Carbon.init(ρC)
  let Pt = Platinum.init(ρPt)
  let glass = compound((Si, 1), (O, 2), ρ = ρGlass) # pure quartz glass

  let m1 = initDepthGradedMultilayer(
    substrate = glass,
    top = C, bottom = Pt,
    dMin = 11.5.nm, dMax = 22.5.nm,
    Γ = 0.45,
    c = 1.0,
    N = 2,
    σ = 0.45.nm
  )
  let m2 = initDepthGradedMultilayer(
    substrate = glass,
    top = C, bottom = Pt,
    dMin = 7.0.nm, dMax = 19.nm,
    Γ = 0.45,
    c = 1.0,
    N = 3,
    σ = 0.45.nm
  )
  let m3 = initDepthGradedMultilayer(
    substrate = glass,
    top = C, bottom = Pt,
    dMin = 5.5.nm, dMax = 16.nm,
    Γ = 0.4,
    c = 1.0,
    N = 4,
    σ = 0.45.nm
  )
  let m4 = initDepthGradedMultilayer(
    substrate = glass,
    top = C, bottom = Pt,
    dMin = 5.0.nm, dMax = 14.0.nm,
    Γ = 0.4,
    c = 1.0,
    N = 5,
    σ = 0.45.nm
  )
  @[m1, m2, m3, m4]

import nimhdf5
from std/os import `/`
proc automatic() =
  let ms = getLayers()
  const θ = 0.5.°
  echo ms[0].reflectivity(0.5.°, 9.0.keV, parallel = false)
  #if true: quit()

  let Energy = 8.047.keV
  let θs = linspace(0.0, 3.0, 1000).mapIt(it.°)

  for i, m in ms:
    var RMs = newSeq[float]()
    for i, θl in θs:
      # 1. compute the refractive indices for each layer
      RMs.add m.reflectivity(θl, Energy, parallel = true)

    ggplot(toDf({"θ" : θs.mapIt(it.float), "R" : RMs}), aes("θ", "R")) +
      geom_line() +
      scale_y_log10() +
      ggtitle(&"Reflectivity of LLNL recipe {i} against angle") +
      ggsave(&"/tmp/llnl_angle_scan_{i}.pdf")

  let energies = linspace(0.03, 10.0, 1000).mapIt(it.keV)
  for i, m in ms:
    var RMs = newSeq[float]()
    for i, E in energies:
      # 1. compute the refractive indices for each layer
      RMs.add m.reflectivity(θ, E, parallel = true)

    #if true: quit()
    ggplot(toDf({"E" : energies.mapIt(it.float), "R" : RMs}), aes("E", "R")) +
      geom_line() +
      ggtitle(&"Reflectivity of LLNL recipe {i}") +
      ggsave(&"/tmp/llnl_reflectivity_auto_{i}.pdf")


  echo ms[0].reflectivity(1.53.°, Energy, parallel = true)
  echo ms[0].reflectivity(1.3.°, Energy, parallel = true)
  #if true: quit()



  #if true: quit()
proc writeH5File() =
  const energyMin = 0.03.keV
  const energyMax = 15.0.keV
  const numEnergy = 1000
  const angleMin = 0.0.°
  const angleMax = 15.0.°
  const numAngle = 1000
  const outfile = "llnl_layer_reflectivities.h5"
  const outpath = "/tmp/"
  let Energy = linspace(energyMin.float, energyMax.float, numEnergy) # scan in a range wider than we use in the raytracer
  let θs = linspace(angleMin.float, angleMax.float, numAngle) # [0.5] # the angle under which we check

  var h5f = H5open(outpath / outfile, "rw")
  var eDset = h5f.create_dataset("Energy", (Energy.len, 1), dtype = float64)
  var aDset = h5f.create_dataset("Angles", (θs.len, 1), dtype = float64)
  eDset[eDset.all] = Energy
  aDset[aDset.all] = θs

  let ms = getLayers()
  for i, m in ms:
    echo "Starting m ", i
    var df = newDataFrame()
    var refl = zeros[float]([θs.len, Energy.len])
    for i, θ in θs:
      for j, E in Energy:
        # 1. compute the refractive indices for each layer
        refl[i, j] = m.reflectivity(θ.°, E.keV, parallel = false)
      let dfLoc = seqsToDf({"Energy" : Energy, "Theta" : θ, "Ref" : refl[i, _].squeeze()})
      df.add dfLoc
    ggplot(df, aes("Energy", "Theta", fill = "Ref")) +
      geom_raster() +
      ggtitle("Reflectivity of LLNL layer type " & $i) +
      ggsave("/tmp/layer" & $i & ".pdf")
    var rDset = h5f.create_dataset("Reflectivity" & $i, (θs.len, Energy.len),
                                   dtype = float64)
    rDset.unsafeWrite(refl.dataArray, refl.size)

  discard h5f.close()

when isMainModule:
  # initalize the elements of the substrate
  manual()
  automatic()
  writeH5File()
