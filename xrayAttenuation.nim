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
    name*: string # name of the element
    nProtons*: int ## runtime value of `Z`. Different name to not clash with `Z` static
    nistDf*: DataFrame # stores lines, μ/ρ (raw data from NIST TSV file)
    nistFormFactorDf*: DataFrame # stores lines, μ/ρ (raw data from NIST TSV file)
    henkeDf*: DataFrame # stores f1, f2 form factors (from Henke TSV files)
    μInterp*: InterpolatorType[float]
    f1Henke*: InterpolatorType[float]
    f2Henke*: InterpolatorType[float]
    f1Nist*: InterpolatorType[float]
    f2Nist*: InterpolatorType[float]
    molarMass*: g•mol⁻¹
    chemSym*: string # the chemical symbol, i.e. shortened name. Read from cache table
    ρ*: g•cm⁻³

  Element*[Z: static int] = ElementRT

  ## NOTE: currently we don't support "recursive compound notation", i.e. (CH₃)₃CH
  ## for the functional groups. Instead use the expanded version!
  NumberElement = tuple[e: ElementRT, number: int]

  Compound* = object
    commonName: string # common name of the compound ("Water", "Salt", ...). Has to be user given
    elements: seq[NumberElement]
    ρ*: g•cm⁻³

  GasMixture* = object
    temperature*: Kelvin
    gases*: seq[Compound]
    ratios*: seq[float] # percentage wise ratios of each gas, i.e. partial pressures from ratio * pressure
    pressure*: MilliBar

  FluorescenceLine* = object
    name*: string
    energy*: keV
    intensity*: float # relative intesity compared to *other lines of the same shell*

iterator pairs*(c: Compound): (ElementRT, int) =
  for (el, num) in c.elements:
    yield (el, num)

iterator pairs*(gm: GasMixture): (Compound, float) =
  for (c, ratio) in zip(gm.gases, gm.ratios):
    yield (c, ratio)

const Resources = currentSourcePath().parentDir() / "xrayAttenuation" / "resources/"
const NIST = "nist_mass_attenuation"
const Henke = "henke_form_factors"
const NIST_scattering_factors = "nist_form_factors"
const ElementTable = CacheTable"Elements"
const ElementSymbolTable = CacheTable"ElementSymbols"
const ElementChemToNameTable = CacheTable"ElementChemToName"
const ElementSeq = CacheSeq"ElementSeq"

when defined(staticBuild):
  ## If this is defined, we will slurp all files in the ~resources~ directory and use those
  ## instead of reading them at RT!
  import std / [os, tables, pathnorm]
  export tables.`[]`, normalizePath
  import shell
  proc buildDataTable(path: string): Table[string, string] =
    ## Reads data of all files in `path`, stores them as (path `string` -> data `string`) mapping
    ##
    ## NOTE: When building with `-d:mingw` the `walkDirRec` is broken. Thus we use `shell`.
    result = initTable[string, string]()
    let Path = Resources.normalizePath(dirSep = '/') #
    let arg = "-name \"*.*\"" # only files!
    let (res, err) = shellVerbose:
      find ($Path) ($arg)
    let files = res.splitLines()
    for x in files:
      result[normalizePath(x, dirSep = '/')] = staticRead(x)
  const dataTable* = buildDataTable(Resources)

macro readDataFile(args: varargs[untyped]): untyped =
  ## Dispatches to either reading from the CT parsed data strings stored in `DataTable`
  ## or a call to `readCsv` otherwise.
  if defined(staticBuild):
    result = nnkCall.newTree(ident"parseCsvString")
    let path = args[0]
    let lookup = nnkBracketExpr.newTree(
      ident"dataTable",
      nnkCall.newTree(ident"normalizePath", path, newLit '/')
    )
    result.add lookup
    for i in 1 ..< args.len:
      result.add args[i]
  else:
    result = nnkCall.newTree(ident"readCsv")
    for arg in args:
      result.add arg
  #echo result.repr

## NOTE: instead of reading this into a global variable, we could also have an option to
## only read it when necessary (to look up compound densities)
let CompoundDensityDf = readDataFile(Resources / "density_common_materials.csv",
                                     sep = ',', header = "#")
let XrayFluroscenceDf = readDataFile(Resources / "xray_line_intensities.csv")
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
  result.add nnkTypeDef.newTree(nnkPostfix.newTree(ident"*", ident"AnyElement"),
                                newEmptyNode(),
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
  AnyCompound* = Compound | AnyElement | ElementRT

type
  ## XXX: If we try to restrict `T` and `B` to `AnyCompound` then we get a
  ## "cannot instantiate" error...
  DepthGradedMultilayer*[T; B; S: AnyCompound] = object
    N*: int ## Number of layer repetitions
    c*: float ## Power used in calculation of thicknesses
    layers*: seq[NanoMeter] ## Thicknesses of each layer
    dMin*: NanoMeter ## Thinnest multilayer (bottom)
    dMax*: NanoMeter ## Thickest multilayer (top)
    Γ: float ## Ratio of the top to bottom layer
    top*: T ## Material at the top of each multilayer
    bottom*: B ## Material at the bottom of each multilayer
    substrate*: S ## Substrate below all N multilayers
    σ*: NanoMeter ## Surface roughness in NanoMeter

proc name*(e: AnyElement | typedesc[AnyElement]): string
proc Z*(e: AnyElement): int

proc readNistData(element: string, z: int): DataFrame =
  let name = if element == "Aluminium": "Aluminum" else: element
  let path = Resources / NIST / &"data_element_{name}_Z_{z}.csv"
  result = readDataFile(path, sep = '\t')
    .mutate(f{float: "Energy[keV]" ~ idx("Energy[MeV]").MeV.to(keV).float})

proc readNistData(element: var AnyElement) =
  element.nistDf = readNistData(element.name(), Z(element))

proc readNistFormFactorData(chemSym: string): DataFrame =
  let path = Resources / NIST_scattering_factors / &"data_element_{chemSym}.csv"
  result = readDataFile(path, sep = ',')

proc readNistFormFactorData(element: var AnyElement) =
  element.nistFormFactorDf = readNistFormFactorData(element.chemSym)

proc readHenkeData(chemSym: string): DataFrame =
  let pathHenke = Resources / Henke / chemSym.toLowerAscii() & ".nff"
  try: # first try with spaces (almost all files)
    result = readDataFile(pathHenke, sep = ' ', skipLines = 1, colNames = @["Energy", "f1", "f2"])
      .rename(f{"Energy[eV]" <- "Energy"})
      .mutate(f{float: "Energy[keV]" ~ idx("Energy[eV]").eV.to(keV).float})
  except IOError:
    # try with tab
    result = readDataFile(pathHenke, sep = '\t', skipLines = 1, colNames = @["Energy", "f1", "f2"])
      .rename(f{"Energy[eV]" <- "Energy"})
      .mutate(f{float: "Energy[keV]" ~ idx("Energy[eV]").eV.to(keV).float})

proc readHenkeData(element: var AnyElement) =
  element.henkeDf = readHenkeData(element.chemSym)

proc readMolarMasses*(): DataFrame =
  result = readDataFile(Resources / "molar_masses.csv", sep = ' ')
    .head(112) # drop last 2 elements
    .mutate(f{Value -> float: "AtomicWeight[g/mol]" ~ (
      if idx("AtomicWeight[g/mol]").kind == VString:
        idx("AtomicWeight[g/mol]").toStr()[1 ..< ^1].parseFloat
      else:
        idx("AtomicWeight[g/mol]").toFloat)
    })

proc f1eval*(it: AnyElement, val: keV): float =
  if val < 0.03.keV:
    result = 0.0
  elif val < 30.keV: # NIST coarse data starts at 2 keV, but we use the range of Henke data
    result = it.f1Henke.eval(val.float)
  else:
    result = it.f1Nist.eval(val.float)

proc f2eval*(it: AnyElement, val: keV): float =
  if val < 0.03.keV:
    result = 0.0
  elif val < 30.keV: # NIST coarse data starts at 2 keV, but we use the range of Henke data
    result = it.f2Henke.eval(val.float)
  else:
    result = it.f2Nist.eval(val.float)

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
  elif result.startsWith("ElementRT"):
    result = $e.chemSym

proc Z*(e: AnyElement): int =
  ## Return the proton number of the given element. This is just the static
  ## int that defines the type.
  e.Z

macro withElements(name, varIdent, body: untyped): untyped =
  ## Constructs a `case` statement that dispatches the given `name` to a
  ## branch in which the element is available as `varIdent`. Essentially:
  ##
  ## case element
  ## of "Hydrogen", "H":
  ##   let typ = Hydrogen.init()
  ##   `body`
  ## ...
  ## else: doAssert false
  result = nnkCaseStmt.newTree(name)
  for el in ElementSeq:
    let elName = ident(el)
    let bodyIt = quote do: # body plus injected type
      let `varIdent` = `elName`.init()
      `body`
    result.add nnkOfBranch.newTree(newLit(el.strVal), # Name of element
                                   newLit(ElementSymbolTable[el.strVal].strVal), # chemical symbol
                                   bodyIt)
  let elseBody = quote do:
    doAssert false, "The element " & $(`name`) & " is not known."
  result.add nnkElse.newTree(elseBody)

macro protonsImpl(name: untyped): untyped =
  ## Constructs a `case` statement that dispatches to the number of protons
  ## for each branch (corresponding to an element or chemical symbol)
  ##
  ## case element
  ## of "Hydrogen", "H":
  ##   result = 1
  ## ...
  ## else: doAssert false
  result = nnkCaseStmt.newTree(name)
  for (el, num) in pairs(ElementTable):
    let elName = ident(el)
    let resId = ident"result"
    let csym = newLit(ElementSymbolTable[el].strVal)
    let bodyIt = quote do: # body plus injected type
      `resId` = `num`
    result.add nnkOfBranch.newTree(newLit(el), # Name of element
                                   csym, # chemical symbol
                                   bodyIt)
  let elseBody = quote do:
    doAssert false, "The element " & $(`name`) & " is not known."
  result.add nnkElse.newTree(elseBody)

macro chemImpl(name: untyped): untyped =
  ## Constructs a `case` statement that dispatches to the chemical symbol
  ## for each branch (corresponding to an element or chemical symbol)
  ##
  ## case element
  ## of "Hydrogen", "H":
  ##   result = "H"
  ## ...
  ## else: doAssert false
  result = nnkCaseStmt.newTree(name)
  for (el, num) in pairs(ElementTable):
    let elName = ident(el)
    let resId = ident"result"
    let csym = newLit(ElementSymbolTable[el].strVal)
    let bodyIt = quote do: # body plus injected type
      `resId` = `csym`
    result.add nnkOfBranch.newTree(newLit(el), # Name of element
                                   csym, # chemical symbol
                                   bodyIt)
  let elseBody = quote do:
    doAssert false, "The element " & $(`name`) & " is not known."
  result.add nnkElse.newTree(elseBody)

macro nameImpl(name: untyped): untyped =
  ## Constructs a `case` statement that dispatches to the name
  ## for each branch (corresponding to an element or chemical symbol)
  ##
  ## case element
  ## of "Hydrogen", "H":
  ##   result = "Hydrogen"
  ## ...
  ## else: doAssert false
  result = nnkCaseStmt.newTree(name)
  for (el, num) in pairs(ElementTable):
    let elName = ident(el)
    let resId = ident"result"
    let csym = newLit(ElementSymbolTable[el].strVal)
    let bodyIt = quote do: # body plus injected type
      `resId` = `el`
    result.add nnkOfBranch.newTree(newLit(el), # name of element
                                   csym, # chemical symbol
                                   bodyIt)
  let elseBody = quote do:
    doAssert false, "The element " & $(`name`) & " is not known."
  result.add nnkElse.newTree(elseBody)
  #echo result.repr

proc protons*(element: string): int =
  ## Returns the number of protons for a given element name or chemical formula
  protonsImpl(element)

proc chemicalSymbol*(element: string): string =
  ## Returns the chemical symbol of the given element
  chemImpl(element)

proc elementName*(element: string): string =
  ## Returns the name of the given element (which may be a chemical symbol!)
  nameImpl(element)

proc initElement*(element: string, ρ = -1.g•cm⁻³): ElementRT =
  ## Return an instance of the desired element, which means reading the
  ## data from the `resources` directory and creating the interpolator to
  ## interpolate arbitrary energies.
  let name = elementName(element)
  result.name = name
  result.nProtons = element.protons()
  result.ρ = ρ
  result.chemSym = chemicalSymbol(element)

  # Read data
  result.nistDf = readNistData(name, result.nProtons)
  result.nistFormFactorDf = readNistFormFactorData(result.chemSym)
  result.henkeDf = readHenkeData(result.chemSym)
  # fill molar mass from global `MolarMassDf`
  result.molarMass = MolarMassDf.filter(f{`Name` == name})["AtomicWeight[g/mol]", float][0].g•mol⁻¹

  result.f1Henke = newLinear1D(result.henkeDf["Energy[keV]", float].toSeq1D,
                               result.henkeDf["f1", float].toSeq1D)
  result.f2Henke = newLinear1D(result.henkeDf["Energy[keV]", float].toSeq1D,
                               result.henkeDf["f2", float].toSeq1D)

  result.f1Nist = newLinear1D(result.nistFormFactorDf["E [keV]", float].toSeq1D,
                              result.nistFormFactorDf["f1 [e atom⁻¹]", float].toSeq1D)
  result.f2Nist = newLinear1D(result.nistFormFactorDf["E [keV]", float].toSeq1D,
                              result.nistFormFactorDf["f2 [e atom⁻¹]", float].toSeq1D)

  ## fix interp!
  ## NIST data has the same energy at every value right _before_ and _on_ a
  ## transition line. Instead we need to modify it such that the value before
  ## is at a slightly lower energy!
  ## So create a lagged column and modify the energy column to have a small shift in
  ## each point that is duplicated
  result.nistDf["EnergyLag"] = lag(result.nistDf["Energy[keV]"])
  result.nistDf = result.nistDf
    .mutate(f{float: "Energy[keV]" ~ (
      if idx("EnergyLag") == idx("Energy[keV]"):
        idx("Energy[keV]") + 0.01 * idx("Energy[keV]")
      else:
        idx("Energy[keV]")
        )})
  result.μInterp = newLinear1D(result.nistDf["Energy[keV]", float].toSeq1D,
                               result.nistDf["μ/ρ", float].toSeq1D)

proc init*[T: AnyElement](element: typedesc[T], ρ = -1.g•cm⁻³): T =
  ## Return an instance of the desired element, which means reading the
  ## data from the `resources` directory and creating the interpolator to
  ## interpolate arbitrary energies.
  let res = initElement($element, ρ)
  result = T(res)

proc name*(c: Compound): string
proc initCompound*(ρ: g•cm⁻³, elements: varargs[(ElementRT, int)]): Compound =
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
  if ρ > 0.g•cm⁻³:
    result.ρ = ρ
  else:
    let densityDf = CompoundDensityDf.filter(f{`Formula` == compoundName})
    if densityDf.len > 0:
      ## NOTE: there Quartz and Silica have the same formula, but different densities.
      ## Currently we take the *first* density!
      result.ρ = densityDf["Density[g•cm⁻³]", float][0].g•cm⁻³

proc parseCompound*(s: string): seq[(ElementRT, int)] =
  ## Parses a given compound string, i.e. `CO2`, `H2O`, `Si3N4`, ...
  var i = 0
  var chemSym = ""
  var num = ""
  var inDigits = false
  while i < s.len:
    case s[i]
    of {'A' .. 'Z'}:
      if inDigits: # last was digit, add last element / number pair
        result.add (initElement(chemSym), parseInt(num))
        chemSym = ""
        num = ""
        inDigits = false
      elif chemSym.len > 0 and num.len == 0: ## Short form, e.g. `CO2` with implied `1`
        result.add (initElement(chemSym), 1)
        chemSym = ""

      chemSym.add s[i]
    of {'a' .. 'z'}: # element with multiple letters
      chemSym.add s[i]
    of {'0' .. '9'}:
      num.add s[i]
      inDigits = true
    else:
      doAssert false, "Encountered unexpected character in compound formula: `" & $s[i] & "`, full input: " & $s
    inc i
  # add last
  if num.len > 0:
    result.add (initElement(chemSym), parseInt(num))
  else:
    result.add (initElement(chemSym), 1)

proc initCompound*(name: string): Compound =
  ## Initializes a compound from a string, i.e. `H2O`, `CO2`, `Si3N4`, ...
  result = parseCompound(name)

proc initGasMixture*[P: Pressure](
  T: Kelvin,
  pressure: P,
  gases: seq[AnyCompound],
  ratios: seq[float]
                                ): GasMixture =
  var sum = 0.0
  for r in ratios:
    sum += r
  if abs(1.0 - sum) > 1e-3:
    raise newException(ValueError, "Given gas mixture does not sum to 100% for " &
      "all contributions: " & $gases)
  result = GasMixture(temperature: T,
                      pressure: pressure.to(MilliBar),
                      ratios: ratios,
                      gases: gases)

proc initGasMixture*[P: Pressure](
  T: Kelvin,
  pressure: P,
  gases: varargs[(AnyCompound, float)]): GasMixture =
  var gs = newSeq[Compound]()
  var ps = newSeq[float]()
  for (name, p) in gases:
    gs.add name
    ps.add p
  result = initGasMixture(T, pressure, gs, ps)

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
  ##
  ## In addition a density can be given
  ##
  ## .. code-block:: Nim
  ##  let H₂O = initCompound((H, 2), (O, 1), ρ = 1.g•cm⁻³)
  ##  # or as
  ##  let H₂O = initCompound((H, 2), (O, 1), density = 1.g•cm⁻³)
  ##
  ## the latter is useful in cases where the compound is not a well known
  ## and tabulated value (in our CSV file).
  var variables = newStmtList()
  var compoundCall = nnkCall.newTree(ident"initCompound")
  # 0. Add the default argument for the density
  let arg = quote do:
    -1.g•cm⁻³
  compoundCall.add nnkExprEqExpr.newTree(ident"ρ", arg)
  var ρ: NimNode
  for arg in args:
    var number = 1 # if no tuple with number of atoms given, default to 1
    var chemSym: string
    case arg.kind
    of nnkTupleConstr:
      doAssert arg[0].kind in {nnkIdent, nnkSym}
      doAssert arg[1].kind == nnkIntLit
      chemSym = arg[0].strVal
      number = arg[1].intVal.int
    of nnkExprEqExpr:
      let argStr = getArgStr(arg[0])
      const allowedArgs = ["ρ", "density"]
      if argStr notin allowedArgs:
        error("Invalid argument: " & $argStr & "!")
      else:
        ρ = arg[1]
      # skip the rest of this iteration
      continue
    of nnkIdent:
      chemSym = arg.strVal
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
  # 3. overwrite the density argument if given
  if ρ.kind != nnkNilLit:
    compoundCall[1] = ρ
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

proc pretty*(gm: GasMixture, showPrefix: bool): string =
  if showPrefix:
    result = "GasMixture: "
  for i, x in gm.gases:
    result.add &"{x} ({gm.ratios[i]})"
    if i < gm.gases.len:
      result.add ", "
  result.add &"(T = {gm.temperature}, P = {gm.pressure})"
proc `$`*(gm: GasMixture): string = pretty(gm, true)

proc initDepthGradedMultilayer*[L: Length; T: AnyCompound; B: AnyCompound; S: AnyCompound](
  top: T, bottom: B, substrate: S,
  dMin, dMax: L,
  Γ: float, N: int, c: float,
  σ: L): DepthGradedMultilayer[T, B, S] =
  let layersMulti = depthGradedLayers(dMin.to(nm), dMax.to(nm), N, c)
  ## Compute the actual layers from Γ and `layers`
  var layers = newSeq[NanoMeter]()
  for layer in layersMulti:
    layers.add layer * Γ
    layers.add layer * (1.0 - Γ)
  result = DepthGradedMultilayer[T, B, S](
    N: N, c: c, Γ: Γ,
    dMin: dMin.to(nm), dMax: dMax.to(nm),
    layers: layers,
    substrate: substrate,
    top: top,
    bottom: bottom,
    σ: σ.to(nm)
  )

proc name*(c: Compound): string = $c

################################
## Physics related procedures ##
################################
proc molarMass*(c: Compound): g•mol⁻¹ =
  ## Computes the molar mass of the compound:
  ##
  ## `M = Σ_i x_i · A_i`
  ## where `x_i` is the number of atoms of that type in the compound
  ## and `A_i` the atomic mass of each atom in `g•mol⁻¹`
  ## (with `A_i = A_r,i · M_u` where `A_r,i` is the standard atomic weight
  ## of the element and `M_u =~= 1 g•mol⁻¹` the molar mass constant).
  for el, num in c:
    result += el.molarMass * num.float

proc molarWeight*(c: Compound): g•mol⁻¹ {.deprecated: "Please use `molarMass(Compound)` instead.".} =
  molarMass(c)

proc molarMass*(gm: GasMixture): g•mol⁻¹ =
  ## Computes the molar weight of a gas mixture.
  for c, ratio in gm:
    result += c.molarWeight() * ratio

proc numAtoms*(c: Compound): int =
  ## Returns the number of atoms in the compound.
  for _, num in c:
    inc result, num

proc atomicAbsorptionCrossSection*(el: AnyElement, energy: keV): cm² =
  ## Computes the atomic absoprtion cross section `σ_a` based on the scattering factor `f2`
  ## via
  ##  `σ_A = 2 r_e λ f₂`
  ## for the given element.
  let λ = wavelength(energy)
  result = (2 * r_e * λ * el.f2eval(energy)).to(cm²)

proc attenuationCoefficient*(e: AnyElement, energy: keV): cm²•g⁻¹ =
  ## Computes the mass attenuation coefficient `μ` for the given element at the
  ## given `energy`
  if energy <= 30.keV:
    result = attenuationCoefficient(energy, e.f2eval(energy), e.molarMass)
  else:
    result = e.μInterp.eval(energy.float).cm²•g⁻¹

proc attenuationCoefficient*(c: Compound, energy: keV): cm²•g⁻¹ =
  ## Computes the mass attenuation coefficient of a `Compound c` at given `energy`
  var factor = N_A / molarWeight(c)
  var sum_σs: cm²
  for el, num in c:
    sum_σs += num.float * el.atomicAbsorptionCrossSection(energy)
  result = factor * sum_σs

proc attenuationCoefficient*(gm: GasMixture, energy: keV): cm⁻¹ =
  ## Computes the attenuation coefficient (multiplied with the density!)
  ## for the gas mixture taking into account the partial pressures of each
  ## gas.
  for (g, r) in zip(gm.gases, gm.ratios):
    let pr = gm.pressure * r # partial pressure of this gas
    let ρr = density(pr, gm.temperature, g.molarWeight())
    result += attenuationCoefficient(g, energy) * ρr

proc transmission*[L: Length, D: Density](c: AnyCompound, ρ: D, length: L, E: keV): float =
  ## Computes the transmission using the Beer-Lambert law of photons of energy `E`
  ## through the given compound with density `ρ` and `length`.
  let μ = c.attenuationCoefficient(E) # attenuation coefficient of the compound
  result = transmission(μ, ρ.to(g•cm⁻³), length.to(Meter))
  if classify(result) == fcNaN and E < 0.03.keV:
    result = c.transmission(ρ, length, 0.03.keV)

proc transmission*[L: Length](c: AnyCompound, length: L, E: keV): float =
  ## Computes the transmission using the Beer-Lambert law of photons of energy `E`
  ## through the density `ρ` of the compound and given `length`.
  if c.ρ == 0.0.g•cm⁻³:
    raise newException(ValueError, "The given compound : " & $c & " does not have a density " &
      "assigned.")
  result = c.transmission(c.ρ, length, E)

proc density*(gm: GasMixture): g•cm⁻³ =
  ## Returns the density of the given gas mixture
  for (g, r) in zip(gm.gases, gm.ratios):
    let pr = gm.pressure * r # partial pressure of this gas
    let ρr = density(pr, gm.temperature, g.molarWeight())
    result += ρr

template ρ*(gm: GasMixture): g•cm⁻³ = gm.density()

proc transmission*[L: Length](gm: GasMixture, length: L, E: keV): float =
  ## Computes the transmission using the Beer-Lambert law of photons of energy `E`
  ## through the given gas mixture with density `ρ` and `length`.
  let μ = gm.attenuationCoefficient(E) # attenuation coefficient of the compound
  let ρ = density(gm)
  ## XXX: Verify that this is actually correct! That we can just use the combined
  ## attenuation coefficient like this! Looks correct for Ar/Iso, but that may be
  ## due to the low fraction on isobutane!
  result = transmission(μ / ρ, ρ.to(g•cm⁻³), length.to(Meter))
  if classify(result) == fcNaN and E < 0.03.keV:
    result = gm.transmission(length, 0.03.keV)

proc absorptionLength*(c: AnyCompound, ρ: g•cm⁻³, energy: keV): Meter =
  ## Computes the absorption length of the given compound and density at `energy`.
  ##
  ## Equivalent to the inverse attenuation cofficient.
  result = (1.0 / (attenuationCoefficient(c, energy) * ρ)).to(Meter)

proc absorptionLength*(gm: GasMixture, energy: keV): Meter =
  ## Computes the absorption length of the given gas mixture at `energy`.
  ##
  ## Equivalent to the inverse attenuation cofficient.
  result = (1.0 / (attenuationCoefficient(gm, energy))).to(Meter)

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

proc refractiveIndex*(c: Compound, energy: keV, ρ: g•cm⁻³): Complex[float] =
  ## Computes the refractive index for the given element at `energy` and density `ρ`.
  ##
  ##  `n(ω) = 1 - r_e λ² / 2π Σ_i n_ai fi(0)`
  ##  `     = 1 - (β + iδ)`
  ##
  ## Computes the second term for each species in the compound first, then
  ## returns `n`.
  ##
  ## We compute it by adding the `(β + iδ)` terms for each element in the compound
  ## weighted by the fractional number density, i.e.
  ## `n_i = n_a * n_atoms * x_i * m_i / M`
  ## where `n_a` is the total number density in *molecules* for this compound,
  ## `n_atoms` the number of atoms in the compound, `x_i` the number of atoms
  ## for the element `i`, `m_i` the molar mass of that element and `M` the total
  ## molar weight of the compound (sum of all molar masses * number of each atoms).
  ##
  ## Important note: This refractive index is of course only valid in the X-ray
  ## regime, because there we do not have to consider the properties of the
  ## molecular and atomic interactions having an effect on the much longer wavelength
  ## in case of visible light.
  doAssert ρ > 0.g•cm⁻³, "Input " & $c & " does not have a density!"
  var sumβδ = complex(0.0, 0.0)
  let numAtoms = numAtoms(c) # total number of atoms in the molecule
  for el, num in c:
    let f0 = el.f0eval(energy) # f0 of this element
    let n_a = numberDensity(ρ, molarWeight(c)) # number density of the molecule
    let fraction = el.molarMass * num / molarWeight(c) # fraction of particles of this type
    let n_i = n_a * numAtoms * fraction
    let n = refractiveIndex(energy, n_i, f0) # refractive index contribution of this element
    sumβδ += (1.0 - n) # correct `1 - (β + iδ)` computed in `refractiveIndex` to get `β + iδ`.
  result = 1.0 - sumβδ # compute back `n` from `1 - (β + iδ)`.

proc reflectivity*(e: AnyElement, energy: keV, ρ: g•cm⁻³, θ: Degree, σ: Meter,
                   parallel: bool): float =
  ## Computes the reflectivity of the given element `e` at the boundary of vacuum to
  ## a flat surface of `e` at the given `energy` and density `ρ`.
  ##
  ## Computed for `p` polarization if `parallel = true`, else `s`-pol.
  ##
  ## `σ` is the surface roughness and is the deviation from a perfectly smooth surface,
  ## approximated by use of the Névot–Croce factor:
  ##  `exp(-2 k_iz k_jz σ²)`
  ## where `k_iz`, `k_jz` are the wave vectors perpendicular to the surface in the medium
  ## before and after the interface between them.
  doAssert ρ > 0.g•cm⁻³, "Input " & $e & " does not have a density!"
  let n = e.refractiveIndex(energy, ρ)
  result = reflectivity(θ, energy, n, σ, parallel = parallel)

proc reflectivity*(e: AnyElement, energy: keV, ρ: g•cm⁻³, θ: Degree, σ: Meter): float =
  ## Overload of the above, which computes it for unpolarized light.
  doAssert ρ > 0.g•cm⁻³, "Input " & $e & " does not have a density!"
  let n = e.refractiveIndex(energy, ρ)
  result = reflectivity(θ, energy, n, σ)

proc reflectivity*[T; B; S](ml: DepthGradedMultilayer[T, B, S],
                            θ_i: Degree, energy: keV, parallel: bool): float =
  # 1. compute the refractive indices at this energy for each materials layer
  const nVacuum = complex(1.0, 0.0)
  let nTop = refractiveIndex(ml.top, energy, ml.top.ρ)
  let nBot = refractiveIndex(ml.bottom, energy, ml.bottom.ρ)
  let nSub = refractiveIndex(ml.substrate, energy, ml.substrate.ρ)
  var ns = @[nVacuum]
  for _ in 0 ..< ml.layers.len div 2:
    ns.add nTop; ns.add nBot
  ns.add nSub
  # 2. compute the reflectivity
  result = multilayerReflectivity(θ_i, energy, ns, ml.layers, parallel)

proc getFluorescenceLines*[E: AnyElement](e: E): seq[FluorescenceLine] =
  ## Returns the relative intensities of all the fluorescence lines. The intensities
  ## are given in numbers ``relative to other lines of the same shell``.
  ##
  ## E.g. Kα1, Kα2, Kβ1, Kβ2, Kγ1 might have intensities with a maximum of 100,
  ## but then Lα1, Lα2, ... will again use intensities from 100 (or some other
  ## value). So normalizations of the relative intensities of different lines
  ## must be made between the same shell (K, L, M, ...). While the numbers _should_
  ## have a maximum of 100 for the most intense line, this is *not actually the case*!
  ## The data is from the X-ray data booklet:
  ##
  ## https://xdb.lbl.gov/Section1/Table_1-3.pdf
  ##
  ## The relative intensity between K and L lines is typically on the order of
  ## 10 to 1.
  ## See for example:
  ##
  ## https://xdb.lbl.gov/Section1/Sec_1-3.html
  let dfZ = XrayFluroscenceDf.filter(f{`Z` == e.Z})
  for r in dfZ:
    result.add FluorescenceLine(name: r["Line"].toStr,
                                energy: r["Energy [eV]"].toFloat.eV.to(keV),
                                intensity: r["Intensity"].toFloat)


###################################
##### Plotting related procs ######
###################################

proc plotAttenuation*(element: Element,
                      range: (float, float) = (0.0, 1e-2),
                      logLog = false,
                      outpath = "/tmp") =
  ## Create a plot of the attenuation coefficients in `outpath` of the given
  ## element. A log-log plot with both `μ/ρ` and `μ_en/ρ`.
  var df = element.nistDf
    .gather(["μ/ρ", "μ_en/ρ"], "Type", "Value")
  let z = Z(element)
  var plt = ggplot(df, aes("Energy[MeV]", "Value", color = "Type")) +
    geom_line() +
    xlab("Photon energy [MeV]") + ylab("Attenuation coefficient")
  if logLog:
    plt = plt + scale_x_log10() + scale_y_log10()
  else:
    plt = plt + xlim(range[0], range[1])
  plt +
    ggtitle(&"Mass attenuation coefficient for: {element.name()} Z = {z}") +
    ggsave(outpath / &"attenuation_{element.name()}.pdf")

proc plotTransmission*[T: AnyCompound](el: T, ρ: g•cm⁻³, length: Meter,
                                       energyMin = 0.03,
                                       energyMax = 10.0,
                                       outpath = "/tmp") =
  ## Plots the relative transmission of X-rays at different energies for the given
  ## element/compound at the given density.
  let lengthStr = length.to(μm).pretty(precision = 2, short = true)
  let densityStr = ρ.pretty(precision = 3, short = true)
  when T is AnyElement:
    var df = el.henkeDf
      .mutate(f{float: "μ" ~ el.attenuationCoefficient(idx("Energy[keV]").keV).float},
              f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, ρ, length).float},
              f{float: "Abs" ~ 1.0 - `Trans`})
    let z = Z(el)
    let title = &"Transmission for: {el.name()} Z = {z}, length = {lengthStr}, at ρ = {densityStr}"
  else:
    let df = toDf({"Energy[keV]" : linspace(energyMin, energyMax, 1000)})
      .mutate(f{float: "μ" ~ el.attenuationCoefficient(idx("Energy[keV]").keV).float},
              f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, ρ, length).float},
              f{float: "Abs" ~ 1.0 - `Trans`})
    let title = &"Transmission for: {el.name()} length = {lengthStr}, at ρ = {densityStr}"
  ggplot(df, aes("Energy[keV]", "Trans")) +
    geom_line() +
    xlim(energyMin, energyMax) +
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
  var R = newSeq[float](df.len) # unpolarized
  var Rs = newSeq[float](df.len) # s-pol
  var Rp = newSeq[float](df.len) # p-pol
  let E = df["Energy[keV]", float]

  for i in 0 ..< df.len:
    δs[i] = element.delta(E[i].keV, ρ)
    βs[i] = element.beta(E[i].keV, ρ)
    let n = complex((1 - δs[i]), βs[i])
    ns[i] = element.refractiveIndex(E[i].keV, ρ)
    R[i] = element.reflectivity(E[i].keV, ρ, θ, σ)
    Rs[i] = element.reflectivity(E[i].keV, ρ, θ, σ, parallel = false)
    Rp[i] = element.reflectivity(E[i].keV, ρ, θ, σ, parallel = true)
  df["n"] = ns.mapIt(it.abs)
  df["R"] = Rs
  df["Rs"] = Rs
  df["Rp"] = Rp
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

  ggplot(df, aes("Energy[keV]", "R")) +
    geom_line() +
    geom_line(aes = aes(y = "Rs"), color = "yellow") +
    geom_line(aes = aes(y = "Rp"), color = "red") +
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
