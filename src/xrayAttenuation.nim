import std / [strformat, os, macros, macrocache, strutils]
import ggplotnim, numericalnim, unchained


import xrayAttenuation / [physics_utils, macro_utils]
export physics_utils

#[
Note:
We use NIST data files based on mass attenuation coefficients for energies above 30 keV
and Henke files based on form factors `f1, f2` for energies from 10 eV to 30 keV.
]#

type
  Element*[Z: static int] = object
    nistDf*: DataFrame # stores lines, μ/ρ (raw data from NIST TSV file)
    henkeDf*: DataFrame # stores f1, f2 form factors (from Henke TSV files)
    μInterp*: InterpolatorType[float]
    f2*: InterpolatorType[float]
    molarMass*: g•mol⁻¹
    chemSym*: string # the chemical symbol, i.e. shortened name. Read from cache table

const Resources = currentSourcePath().parentDir().parentDir() / "resources"
const NIST = "nist_mass_attenuation"
const Henke = "henke_form_factors"
const ElementTable = CacheTable"Elements"
const ElementSymbolTable = CacheTable"ElementSymbols"
const ElementSeq = CacheSeq"ElementSeq"

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

proc name*(e: AnyElement | typedesc[AnyElement]): string
proc Z*(e: AnyElement): int

proc readNistData(element: var AnyElement) =
  let z = Z(element)
  echo Resources / NIST / &"data_element_{element.name()}_Z_{z}.csv"
  element.nistDf = readCsv(Resources / NIST / &"data_element_{element.name()}_Z_{z}.csv",
                           sep = '\t')
    .mutate(f{float: "Energy[keV]" ~ idx("Energy[MeV]").MeV.to(keV).float})

proc readHenkeData(element: var AnyElement) =
  let path = Resources / Henke / element.chemSym.toLowerAscii() & ".nff"
  element.henkeDf = readCsv(path, sep = ' ', skipLines = 1, colNames = @["Energy", "f1", "f2"])
    .rename(f{"Energy[eV]" <- "Energy"})
    .mutate(f{float: "Energy[keV]" ~ idx("Energy[eV]").eV.to(keV).float})

proc readMolarMasses(): DataFrame =
  result = readCsv(Resources / "molar_masses.csv", sep = ' ')
    .head(112) # drop last 2 elements
    .mutate(f{Value -> float: "AtomicWeight[g/mol]" ~ (
      if idx("AtomicWeight[g/mol]").kind == VString:
        idx("AtomicWeight[g/mol]").toStr()[1 ..< ^1].parseFloat
      else:
        idx("AtomicWeight[g/mol]").toFloat)
    })

proc f2eval*(it: AnyElement, val: keV): float =
  it.f2.eval(val.float)

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
  result.chemSym = lookupChemSymbol(element)
  result.readHenkeData()
  # fill molar mass
  let mM = readMolarMasses()
  result.molarMass = mM.filter(f{`Name` == name})["AtomicWeight[g/mol]", float][0].g•mol⁻¹

  result.f2 = newLinear1D(result.henkeDf["Energy[keV]", float].toSeq1D,
                          result.henkeDf["f2", float].toSeq1D)
  ## fix interp!
  ## NOTE: In order to get it working, we need to modify the mass attenuation coefficients
  ## we download from NIST. It currently has the same energy at every value right _before_
  ## and _on_ a transition line. Instead we need to modify it such that the value before
  ## is at a slightly lower energy!
  #result.μInterp = newLinear1D(result.nistDf["Energy[keV]", float].toSeq1D,
  #                             result.nistDf["μ/ρ", float].toSeq1D)

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

proc plotTransmission*(element: Element, ρ: g•cm⁻³, length: Meter, outpath = "/tmp") =
  ## Plots the relative transmission of X-rays at different energies for the given
  ## element at the given density.
  var df = element.henkeDf
    .mutate(f{float: "μ" ~ attenuationCoefficient(idx("Energy[keV]").keV,
                                                  idx("f2"),
                                                  element.molarMass).float},
            f{float: "Trans" ~ transmission(`μ`.cm²•g⁻¹, ρ, length).float},
            f{float: "Abs" ~ 1.0 - `Trans`})
  let z = Z(element)
  ggplot(df, aes("Energy[keV]", "Trans")) +
    geom_line() +
    xlim(0.0, 10) +
    xlab("Photon energy [keV]") + ylab("Transmission") +
    ggtitle(&"Transmission for: {element.name()} Z = {z}, length = {length}, at ρ = {ρ}") +
    ggsave(outpath / &"transmission_{element.name()}.pdf", width = 800, height = 480)

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

when isMainModule:
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
