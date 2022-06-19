import std / [strformat, os, macros, macrocache, strutils, complex, sequtils]
import ggplotnim, numericalnim, unchained

import xrayAttenuation / [physics_utils, macro_utils]
export physics_utils

# need units in user program
export unchained

#[
Note:
We use NIST data files based on mass attenuation coefficients for energies above 30 keV
and Henke files based on form factors `f1, f2` for energies from 10 eV to 30 keV.
]#

type
  ElementRT* = object
    nProtons*: int ## runtime value of `Z`. Different name to not clash with `Z` static
    nistDf*: DataFrame # stores lines, μ/ρ (raw data from NIST TSV file)
    henkeDf*: DataFrame # stores f1, f2 form factors (from Henke TSV files)
    μInterp*: InterpolatorType[float]
    f1*: InterpolatorType[float]
    f2*: InterpolatorType[float]
    molarMass*: g•mol⁻¹
    chemSym*: string # the chemical symbol, i.e. shortened name. Read from cache table

  Element*[Z: static int] = ElementRT

  ## NOTE: currently we don't support "recursive compound notation", i.e. (CH₃)₃CH
  ## for the functional groups. Instead use the expanded version!
  NumberElement = tuple[e: ElementRT, number: int]

  Compound* = object
    commonName: string # common name of the compound ("Water", "Salt", ...). Has to be user given
    elements: seq[NumberElement]
    ρ: g•cm⁻³


iterator pairs(c: Compound): (ElementRT, int) =
  for (el, num) in c.elements:
    yield (el, num)

const Resources = currentSourcePath().parentDir().parentDir() / "resources"
const NIST = "nist_mass_attenuation"
const Henke = "henke_form_factors"
const ElementTable = CacheTable"Elements"
const ElementSymbolTable = CacheTable"ElementSymbols"
const ElementChemToNameTable = CacheTable"ElementChemToName"
const ElementSeq = CacheSeq"ElementSeq"

## NOTE: instead of reading this into a global variable, we could also have an option to
## only read it when necessary (to look up compound densities)
let CompoundDensityDf = readCsv(Resources / "density_common_materials.csv",
                                sep = ',', header = "#")
proc readMolarMasses*(): DataFrame
let MolarMassDf = readMolarMasses()

macro generateElements(elms: untyped): untyped =
  ## Generate all given elements from `name = Z` pairs and stores the
  ## information in two macrocache instances to lookup information later if needed.
  result = nnkTypeSection.newTree()
  var elements = newSeq[NimNode]()
  for el in elms:
    doAssert el.kind == nnkAsgn
    doAssert el[1].kind == nnkTupleConstr
    let name = el[0]
    let protons = el[1][0]
    let chemSym = el[1][1]
    result.add nnkTypeDef.newTree(nnkPostfix.newTree(ident"*", name),
                                  newEmptyNode(),
                                  nnkBracketExpr.newTree(ident"Element",
                                                         protons))
    ElementTable[name.strVal] = protons
    ElementSymbolTable[name.strVal] = chemSym
    ElementChemToNameTable[chemSym.strVal] = name
    ElementSeq.add name

    elements.add name
  let orN = elements.genTypeClass()
  result.add nnkTypeDef.newTree(ident"AnyElement", newEmptyNode(),
                                orN)

## The following generates types with static int arguments representing their
## atomic number and stores information in CT macrocache objects.
generateElements:
  Hydrogen     = (1,  H )
  Helium       = (2,  He)
  Lithium      = (3,  Li)
  Beryllium    = (4,  Be)
  Boron        = (5,  B )
  Carbon       = (6,  C )
  Nitrogen     = (7,  N )
  Oxygen       = (8,  O )
  Fluorine     = (9,  F )
  Neon         = (10, Ne)
  Sodium       = (11, Na)
  Magnesium    = (12, Mg)
  Aluminium    = (13, Al)
  Silicon      = (14, Si)
  Phosphorus   = (15, P )
  Sulfur       = (16, S )
  Chlorine     = (17, Cl)
  Argon        = (18, Ar)
  Potassium    = (19, K )
  Calcium      = (20, Ca)
  Scandium     = (21, Sc)
  Titanium     = (22, Ti)
  Vanadium     = (23, V )
  Chromium     = (24, Cr)
  Manganese    = (25, Mn)
  Iron         = (26, Fe)
  Cobalt       = (27, Co)
  Nickel       = (28, Ni)
  Copper       = (29, Cu)
  Zinc         = (30, Zn)
  Gallium      = (31, Ga)
  Germanium    = (32, Ge)
  Arsenic      = (33, As)
  Selenium     = (34, Se)
  Bromine      = (35, Br)
  Krypton      = (36, Kr)
  Rubidium     = (37, Rb)
  Strontium    = (38, Sr)
  Yttrium      = (39, Y )
  Zirconium    = (40, Zr)
  Niobium      = (41, Nb)
  Molybdenum   = (42, Mo)
  Technetium   = (43, Tc)
  Ruthenium    = (44, Ru)
  Rhodium      = (45, Rh)
  Palladium    = (46, Pd)
  Silver       = (47, Ag)
  Cadmium      = (48, Cd)
  Indium       = (49, In)
  Tin          = (50, Sn)
  Antimony     = (51, Sb)
  Tellurium    = (52, Te)
  Iodine       = (53, I )
  Xenon        = (54, Xe)
  Caesium      = (55, Cs)
  Barium       = (56, Ba)
  Lanthanum    = (57, La)
  Cerium       = (58, Ce)
  Praseodymium = (59, Pr)
  Neodymium    = (60, Nd)
  Promethium   = (61, Pm)
  Samarium     = (62, Sm)
  Europium     = (63, Eu)
  Gadolinium   = (64, Gd)
  Terbium      = (65, Tb)
  Dysprosium   = (66, Dy)
  Holmium      = (67, Ho)
  Erbium       = (68, Er)
  Thulium      = (69, Tm)
  Ytterbium    = (70, Yb)
  Lutetium     = (71, Lu)
  Hafnium      = (72, Hf)
  Tantalum     = (73, Ta)
  Tungsten     = (74, W )
  Rhenium      = (75, Re)
  Osmium       = (76, Os)
  Iridium      = (77, Ir)
  Platinum     = (78, Pt)
  Gold         = (79, Au)
  Mercury      = (80, Hg)
  Thallium     = (81, Tl)
  Lead         = (82, Pb)
  Bismuth      = (83, Bi)
  Polonium     = (84, Po)
  Astatine     = (85, At)
  Radon        = (86, Rn)
  Francium     = (87, Fr)
  Radium       = (88, Ra)
  Actinium     = (89, Ac)
  Thorium      = (90, Th)
  Protactinium = (91, Pa)
  Uranium      = (92, U )

type
  AnyCompound = Compound | AnyElement

## Usually not a fan of converters, but this should be safe. We only need it to
## construct `Compounds` from any statically typed `Element`
#converter toElementRT*(e: AnyElement): ElementRT = ElementRT(e)

proc name*(e: AnyElement | typedesc[AnyElement]): string
proc Z*(e: AnyElement): int

proc readNistData(element: var AnyElement) =
  let z = Z(element)
  let path = Resources / NIST / &"data_element_{element.name()}_Z_{z}.csv"
  echo path
  element.nistDf = readCsv(path, sep = '\t')
    .mutate(f{float: "Energy[keV]" ~ idx("Energy[MeV]").MeV.to(keV).float})

proc readHenkeData(element: var AnyElement) =
  let path = Resources / Henke / element.chemSym.toLowerAscii() & ".nff"
  echo path
  # first try with spaces (almost all files)
  try:
    element.henkeDf = readCsv(path, sep = ' ', skipLines = 1, colNames = @["Energy", "f1", "f2"])
      .rename(f{"Energy[eV]" <- "Energy"})
      .mutate(f{float: "Energy[keV]" ~ idx("Energy[eV]").eV.to(keV).float})
  except IOError:
    # try with tab
    element.henkeDf = readCsv(path, sep = '\t', skipLines = 1, colNames = @["Energy", "f1", "f2"])
      .rename(f{"Energy[eV]" <- "Energy"})
      .mutate(f{float: "Energy[keV]" ~ idx("Energy[eV]").eV.to(keV).float})

proc readMolarMasses*(): DataFrame =
  result = readCsv(Resources / "molar_masses.csv", sep = ' ')
    .head(112) # drop last 2 elements
    .mutate(f{Value -> float: "AtomicWeight[g/mol]" ~ (
      if idx("AtomicWeight[g/mol]").kind == VString:
        idx("AtomicWeight[g/mol]").toStr()[1 ..< ^1].parseFloat
      else:
        idx("AtomicWeight[g/mol]").toFloat)
    })

proc f1eval*(it: AnyElement, val: keV): float =
  it.f1.eval(val.float)

proc f2eval*(it: AnyElement, val: keV): float =
  it.f2.eval(val.float)

proc f0eval*(it: AnyElement, val: keV): Complex[float] =
  ## Compute the `f0(ω)` value, the scattering factor (forward scattering)
  ##
  ## Note: the sign in the imaginary part depends on the convention used
  ## for the plane wave description `exp(-i(ωt - kr))` vs. `exp(i(ωt - kr))`
  ## where the negative sign corresponds to the latter.
  result = complex(it.f1eval(val), -it.f2eval(val))

proc name*(e: AnyElement | typedesc[AnyElement]): string =
  ## Returns the name of the given element. Takes care of converting
  ## `Element[Z]` style "types" into their correct names.
  result = $typeof(e)
  if result.startsWith("Element["):
    result = $lookupInverseName(e)

proc Z*(e: AnyElement): int =
  ## Return the proton number of the given element. This is just the static
  ## int that defines the type.
  e.Z

proc init*[T: AnyElement](element: typedesc[T]): T =
  ## Return an instance of the desired element, which means reading the
  ## data from the `resources` directory and creating the interpolator to
  ## interpolate arbitrary energies.
  result.readNistData()
  let name = element.name()
  result.nProtons = result.Z
  result.chemSym = lookupChemSymbol(element)
  result.readHenkeData()
  # fill molar mass from global `MolarMassDf`
  result.molarMass = MolarMassDf.filter(f{`Name` == name})["AtomicWeight[g/mol]", float][0].g•mol⁻¹

  result.f1 = newLinear1D(result.henkeDf["Energy[keV]", float].toSeq1D,
                          result.henkeDf["f1", float].toSeq1D)
  result.f2 = newLinear1D(result.henkeDf["Energy[keV]", float].toSeq1D,
                          result.henkeDf["f2", float].toSeq1D)
  ## fix interp!
  ## NOTE: In order to get it working, we need to modify the mass attenuation coefficients
  ## we download from NIST. It currently has the same energy at every value right _before_
  ## and _on_ a transition line. Instead we need to modify it such that the value before
  ## is at a slightly lower energy!
  ## TODO: We don't really need `μInterp` anymore. Just remove?
  #result.μInterp = newLinear1D(result.nistDf["Energy[keV]", float].toSeq1D,
  #                             result.nistDf["μ/ρ", float].toSeq1D)

proc name*(c: Compound): string
proc initCompound*(elements: varargs[(ElementRT, int)]): Compound =
  ## Initializes a `Compound` based on the given `Elements` and the number of
  ## atoms of that kind in the `Compound`.
  ##
  ## Note that the element type given is `ElementRT`, which is the "base"
  ## `Element` type that the specific `Elements` are based on. There is a
  ## `converter` from
  for arg in elements:
    result.elements.add (e: ElementRT(arg[0]), number: arg[1])
  # see if there is a tabulated value for this compound in our common density file
  let compoundName = result.name()
  let densityDf = CompoundDensityDf.filter(f{`Formula` == compoundName})
  if densityDf.len > 0:
    ## NOTE: there Quartz and Silica have the same formula, but different densities.
    ## Currently we take the *first* density!
    result.ρ = densityDf["Density[g•cm⁻³]", float][0].g•cm⁻³

proc initCompound*(name: string): Compound =
  ## TODO: implement a parser for `Formula` in the CompoundDensityDf. Then we can generate a
  ## Compound straight from a name and fill the `name` field!
  doAssert false, "Parsing of `Formula` names not yet implemented!"

macro compound*(args: varargs[untyped]): untyped =
  ## Generates a `Compound` from given the given chemical symbols. If a
  ## tuple is given the first field refers to the element and the second to the
  ## number of atoms of that type.
  ##
  ## The given chemical symbol generates a local instance of the element under
  ## that name and uses it to generate the Compound. Therefore, one does not
  ## need to first `init` an element. Instead of:
  ##
  ## .. code-block:: Nim
  ##  let H = Hydrogin.init()
  ##  let O = Oxygen.init()
  ##  let H₂O = initCompound((H, 2), (O, 1))
  ##
  ## it is enough to write:
  ##
  ## .. code-block:: Nim
  ##  let H₂O = compound (H, 2), O
  var variables = newStmtList()
  var compoundCall = nnkCall.newTree(ident"initCompound")
  for arg in args:
    var number = 1 # if no tuple with number of atoms given, default to 1
    var chemSym: string
    case arg.kind
    of nnkTupleConstr:
      doAssert arg[0].kind == nnkIdent
      doAssert arg[1].kind == nnkIntLit
      chemSym = arg[0].strVal
      number = arg[1].intVal.int
    of nnkIdent:
      let chemSym = arg.strVal
    else:
      error("The `compound` macro receives either a tuple of chemical symbol and number " &
        "of atoms `(H, 2)` or only a chemical symbol `O`.")

    # 1. generate symbols & instanciate the element
    let elName = genSym(nskLet, chemSym)
    let fullName = lookupNameFromChemSymbol(chemSym)
    variables.add quote do:
      let `elName` = `fullName`.init()
    # 2. add the created element and its number of atoms to the `initCompound` call
    compoundCall.add nnkTupleConstr.newTree(elName, newLit number)
  # first init all elements
  result = variables
  # then perform the call
  result.add compoundCall
  # wrap it in a block
  result = quote do:
    block:
      `result`

proc `$`*(c: Compound): string =
  ## Stringifies the compound by generating the chemical formula for it.
  for (el, num) in c.elements:
    if num == 1:
      result.add $el.chemSym
    else:
      result.add $el.chemSym & $num

proc name*(c: Compound): string = $c


################################
## Physics related procedures ##
################################
proc molarWeight*(c: Compound): g•mol⁻¹ =
  ## Computes the molar weight of the compound:
  ##
  ## `MW = Σ_i x_i · A_i`
  ## where `x_i` is the number of atoms of that type in the compound
  ## and `A_i` the atomic mass of each atom in `g•mol⁻¹`.
  for el, num in c:
    result += el.molarMass * num.float

proc atomicAbsorptionCrossSection*(el: AnyElement, energy: keV): cm² =
  ## Computes the atomic absoprtion cross section `σ_a` based on the scattering factor `f2`
  ## via
  ##  `σ_A = 2 r_e λ f₂`
  ## for the given element.
  let λ = wavelength(energy)
  result = (2 * r_e * λ * el.f2eval(energy)).to(cm²)

proc attenuationCoefficient*(e: AnyElement, energy: keV): cm²•g⁻¹ =
  ## Computes the attenuation coefficient `μ` for the given element at the
  ## given `energy`
  result = attenuationCoefficient(energy, e.f2eval(energy), e.molarMass)

proc attenuationCoefficient*(c: Compound, energy: keV): cm²•g⁻¹ =
  ## Computes the attenuation coefficient of a `Compound c` at given `energy`
  var factor = N_A / molarWeight(c)
  var sum_σs: cm²
  for el, num in c:
    sum_σs += num.float * el.atomicAbsorptionCrossSection(energy)
  result = factor * sum_σs

proc delta*(e: AnyElement, energy: keV, ρ: g•cm⁻³): float =
  ## Computes the `delta` of the element at `energy` and density `ρ`
  result = delta(energy, numberDensity(ρ, e.molarMass), e.f1eval(energy))

proc beta*(e: AnyElement, energy: keV, ρ: g•cm⁻³): float =
  ## Computes the `beta` of the element at `energy` and density `ρ`
  result = beta(energy, numberDensity(ρ, e.molarMass), e.f2eval(energy))

proc refractiveIndex*(e: AnyElement, energy: keV, ρ: g•cm⁻³): Complex[float] =
  ## Computes the refractive index for the given element at `energy` and density `ρ`.
  let f0 = e.f0eval(energy)
  result = refractiveIndex(energy, numberDensity(ρ, e.molarMass), f0)

proc reflectivity*(e: AnyElement, energy: keV, ρ: g•cm⁻³, θ: Degree, σ: Meter): float =
  ## Computes the reflectivity of the given element `e` at the boundary of vacuum to
  ## a flat surface of `e` at the given `energy` and density `ρ`. `σ` is the surface
  ## roughness and is the deviation from a perfectly smooth surface, approximated by
  ## use of the Névot–Croce factor:
  ##  `exp(-2 k_iz k_jz σ²)`
  ## where `k_iz`, `k_jz` are the wave vectors perpendicular to the surface in the medium
  ## before and after the interface between them.
  let n = e.refractiveIndex(energy, ρ)
  result = reflectivity(θ, energy, n, σ)

################################
## Plotting related procs ######
################################

proc plotAttenuation*(element: Element, outpath = "/tmp") =
  ## Create a plot of the attenuation coefficients in `outpath` of the given
  ## element. A log-log plot with both `μ/ρ` and `μ_en/ρ`.
  var df = element.nistDf
    .gather(["μ/ρ", "μ_en/ρ"], "Type", "Value")
  let z = Z(element)
  ggplot(df, aes("Energy[MeV]", "Value", color = "Type")) +
    geom_line() +
    xlim(0.0, 1e-2) +
    xlab("Photon energy [MeV]") + ylab("Attenuation coefficient") +
    #scale_x_log10() + scale_y_log10() +
    ggtitle(&"Mass attenuation coefficient for: {element.name()} Z = {z}") +
    ggsave(outpath / &"attenuation_{element.name()}.pdf")

proc plotTransmission*[T: AnyCompound](el: T, ρ: g•cm⁻³, length: Meter, outpath = "/tmp") =
  ## Plots the relative transmission of X-rays at different energies for the given
  ## element/compound at the given density.
  when T is AnyElement:
    var df = el.henkeDf
      .mutate(f{float: "μ" ~ el.attenuationCoefficient(idx("Energy[keV]").keV).float},
              f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, ρ, length).float},
              f{float: "Abs" ~ 1.0 - `Trans`})
    let z = Z(el)
    let title = &"Transmission for: {el.name()} Z = {z}, length = {length}, at ρ = {ρ}"
  else:
    let df = toDf({"Energy[keV]" : linspace(0.03, 10.0, 1000)})
      .mutate(f{float: "μ" ~ el.attenuationCoefficient(idx("Energy[keV]").keV).float},
              f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, ρ, length).float},
              f{float: "Abs" ~ 1.0 - `Trans`})
    let title = &"Transmission for: {el.name()} length = {length}, at ρ = {ρ}"
  ggplot(df, aes("Energy[keV]", "Trans")) +
    geom_line() +
    xlim(0.0, 3.0) +
    xlab("Photon energy [keV]") + ylab("Transmission") +
    ggtitle(title) +
    ggsave(outpath / &"transmission_{el.name()}.pdf", width = 800, height = 480)

proc plotReflectivity*(element: Element, ρ: g•cm⁻³,
                       θ: Degree, σ = 0.0.m,
                       outpath = "/tmp",
                       energyMin = 0.03, energyMax = 10.0): DataFrame =
  ## Plots the reflectivity, refractive index and `δ`, `β` (terms making up refractive index)
  ## of X-rays at different energies for the given element at the given density.
  # remove all data points outside the given range
  var df = element.henkeDf.filter(f{idx("Energy[keV]") > energyMin and idx("Energy[keV]") < energyMax})
  var ns = newSeq[Complex[float]](df.len)
  var δs = newSeq[float](df.len)
  var βs = newSeq[float](df.len)
  var Rs = newSeq[float](df.len)
  let E = df["Energy[keV]", float]
  for i in 0 ..< df.len:
    δs[i] = element.delta(E[i].keV, ρ)
    βs[i] = element.beta(E[i].keV, ρ)
    let n = complex((1 - δs[i]), βs[i])
    ns[i] = element.refractiveIndex(E[i].keV, ρ)
    Rs[i] = element.reflectivity(E[i].keV, ρ, θ, σ)
  df["n"] = ns.mapIt(it.abs)
  df["Rs"] = Rs
  df["δs"] = δs
  df["βs"] = βs

  let z = Z(element)
  ggplot(df, aes("Energy[keV]", "n")) +
    geom_line() +
    xlim(0.0, 10) +
    xlab("Photon energy [keV]") + ylab("Refractive index (absolute value)") +
    ggtitle(&"Refractive index for: {element.name()} Z = {z}, angle = {θ}, at ρ = {ρ}") +
    ggsave(outpath / &"refractive_index_{element.name()}.pdf", width = 800, height = 480)
  showBrowser(df)
  ggplot(df, aes("Energy[keV]")) +
    geom_line(aes = aes(y = "δs", color = "δ")) +
    geom_line(aes = aes(y = "βs", color = "β")) +
    scale_x_log10() + scale_y_log10() +
    xlab("Photon energy [keV]") + ylab("δ, β") +
    ggtitle(&"δ, β for: {element.name()} Z = {z}, angle = {θ}, at ρ = {ρ}") +
    ggsave(outpath / &"delta_beta_index_{element.name()}.pdf", width = 800, height = 480)

  ggplot(df, aes("Energy[keV]", "Rs")) +
    geom_line() +
    #xlim(0.0, 10) +
    xlab("Photon energy [keV]") + ylab("Refractive index (absolute value)") +
    ggtitle(&"Reflectivity for: {element.name()} Z = {z}, angle = {θ}, σ = {σ}, at ρ = {ρ}") +
    ggsave(outpath / &"reflectivity_{element.name()}.pdf", width = 800, height = 480)

  result = df

when false:
  # NOTE: These will become the basis for what is already done in `plotTransmission`
  # if the user wants a specific value instead of a plot.
  proc beerLambert(element: AnyElement,
                   energy: keV,
                   p: Pascal,
                   T: Kelvin,
                   length: CentiMeter) =
    ## Absoprtion for a gas. Density will be computed using ideal gas law.
    ##
    ## The absorption will be computed at the given `energy` and pressure / temperature
    ## conditions.
    ##
    ## This only returns ``relative`` absorption, i.e. assuming `I_0 = 1`.
    let d = density(p, T, element.molarMass)
    #result = exp(-element.

  template absorption(): untyped = beerLambert()
