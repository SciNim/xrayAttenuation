import xrayAttenuation, datamancer, unchained
import std / [strformat, terminal, strutils]
from std / sequtils import mapIt

#[
nim musl -d:release -d:staticBuild calc_density_fraction
]#

proc getMassFractions(c: AnyCompound): seq[(string, float)] =
  ## Computes the mass fractions of all elements inside of a compound.
  for el, num in c:
    let fraction = el.molarMass * num / molarWeight(c) # fraction of particles of this type
    result.add (el.name(), fraction.float)

proc getFractions(gm: GasMixture): seq[(string, float)] =
  ## Computes the particle fractions for a gas mixture (i.e. the fractions used
  ## for partial pressures), *BUT* with the addition that the individual compounds
  ## in a gas, i.e. Isobutane, are split up *BY MASS FRACTION* as well. This is
  ## very confusing and wrong (because ideal gas physics only cares about the
  ## number of *particles*, whether they be atoms or molecules), but this is what
  ## was done in REST in the past.
  for c, ratio in gm:
    for (el, frac) in getMassFractions(c):
      result.add (el, frac * ratio)

proc getMassFractions(gm: GasMixture): seq[(string, float)] =
  ## We compute the gas mass fraction as:
  ##
  ## f_m,k = X_k * M_k / M_gm                  (1)
  ##
  ## with `X_k` the fractional number of particles of that part of the gas
  ## compound, `M_k` the molar mass of the compound in the gas and `M_gm` the
  ## total molar mass of the gas mixture.
  ##
  ## For individual atoms of a compound, i.e. C in Ar 97.7 / Isobutane 2.3 (C4H10)
  ## it is further split by the individual fractions based on the weighted masses
  ## in the compound.
  ##
  ## fm_k = f_m,c * N_c * A_c * M_u / M_c       (2)
  ##
  ## where `f_m,c` is the mass fraction of the compound e.g. Isobutane following eq. (1)
  ## and `N_c` the number of atoms of the compound, e.g. `4` carbon per isobutane
  ## `A_c` the standard atomic weight e.g. `~12.011` for carbon and `M_u` the
  ## molar mass constant (prior 2019 redefinition `1 g/mol`) and `M_c` the
  ## molar mass of the compound (e.g. `58.12 g•mol⁻¹` for Isobutane).
  for c, ratio in gm:
    for (el, cFrac) in getMassFractions(c): # cFrac == compound fraction
      let fraction = cFrac * ratio * molarMass(c) / molarMass(gm)
      result.add (el, float fraction)

proc gas(elements: seq[string], pressure: mbar, T = 293.15.K) =
  var gases: seq[Compound]
  var frac: seq[float]
  for el in elements:
    let sp = el.split(",").mapIt(it.strip)
    let elRTs = parseCompound(sp[0])
    let fr = parseFloat sp[1]

    gases.add initCompound(0.0.g•cm⁻³, elRTs)
    if fr > 1.0:
      echo "[ERROR] Pleas give the gas fractions as relative to 1.0 and not as percentages."
      return
    frac.add fr

  if abs(frac.sum - 1) > 1e-4:
    echo frac.sum
    echo "From numbers: ", frac
    echo "[ERROR] The gas mixture fractions do not sum to 1!"
    return

  let gm = initGasMixture(T, pressure, gases, frac)
  stdout.styledWriteLine(fgYellow, "==================== Gas ====================")
  stdout.styledWriteLine(fgYellow, &"\t{gm}")
  stdout.styledWriteLine(fgYellow, &"\tρ = {gm.ρ.to(kg•m⁻³)}")
  stdout.styledWriteLine(fgRed, "---------- Particle fractions ----------")
  for fr in getFractions(gm):
    stdout.styledWriteLine(fgRed, &"\t{fr[0]:<2}: {fr[1]}")
  stdout.styledWriteLine(fgRed, "------------ Mass fractions ------------")
  for fr in getMassFractions(gm):
    stdout.styledWriteLine(fgRed, &"\t{fr[0]:<2}: {fr[1]}")

proc solid(elements: seq[string], ρ = -1.0.g•cm⁻³) =
  var elRT: seq[(ElementRT, int)]
  for el in elements:
    elRT.add parseCompound(el)

  var comp = initCompound(ρ, elRT)
  stdout.styledWriteLine(fgYellow, "==================== Solid ====================")
  stdout.styledWriteLine(fgYellow, &"\t{comp}")
  stdout.styledWriteLine(fgYellow, &"\tρ = {comp.ρ.to(g•cm⁻³)}")
  stdout.styledWriteLine(fgRed, "------------ Mass fractions ------------")
  for el, fr in getMassFractions(comp):
    stdout.styledWriteLine(fgRed, &"\t{el:<2}: {fr}")

when isMainModule:
  import cligen
  import unchained / cligenParseUnits
  dispatchMulti([gas,
                 help = {
                   "elements" : """Give each part of the gas as a pair 'compound,fraction'
Example: Ar,0.977 C4H10,0.023

For molecular gases (in the example isobutane), simply provide the formula of the molecule.
""",
                   "pressure" : "The pressure of the gas in MilliBar.",
                   "T"        : "The temperature of the gas in Kelvin."
                   }
                ],
                [solid,
                 help = {
                   "elements" : "Give each element part of the compound as a simplified chemical formula. Example: C6H10",
                   "ρ"        : "Density of the compound."}
                 ])
