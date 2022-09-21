* =xrayAttenuation=

This library provides the ability to compute transmission / absorption
behavior of X-rays at different energies for all elements below Z
< 93.

Mass attenuation coefficients μ_m/ρ are used in the Beer-Lambert law
to compute the absorption behavior. For most energies these are first
computed from scattering form factors.

The library basically provides Nim types for all different elements
that are supported via a =Element[Z: static int]= base type and
compounds of multiple elements via a =Compound= type, which stores the
element and number of atoms of that type.

Example:
#+begin_src nim
import xrayAttenuation
  
let ar = Argon.init() # generate an Argon instance
ar.plotAttenuation() # generate a plot of the attenuation factors against energy

# compute a density at known pressure and temperature (the molar mass is filled automatically
# in `init`)
let ρ_Ar = density(1050.mbar.to(Pascal), 293.K, ar.molarMass)
# use our density and a desired length to generate a plot of the transmission in 3cm Argon
ar.plotTransmission(ρ_Ar, 3.cm.to(m))
# compute an individual absorption length at a desired energy
echo absorptionLength(2.5.keV, numberDensity(ρ_Ar, ar.molarMass),
                      ar.f2.eval(2.5))

# this can also be computed for compounds. To make construction of a
# compound easier, we can use the `compound` macro
block SimpleCompound:
  let Si = Silicon.init()
  let N = Nitrogen.init()
  let Si₃N₄ = initCompound((Si, 3), (N, 4))
  Si₃N₄.plotTransmission(Si₃N₄.ρ, # for some common compounds we have a table of densities,
                                 # which are filled upon compound initialization
                         300.nm.to(Meter))
block CompoundMacro:
  # or simpler using the macro (i.e. no initialization of `Si` and `N` needed:
  let Si₃N₄ = compound (Si, 3), (N, 4)
  Si₃N₄.plotTransmission(Si₃N₄.ρ, 300.nm.to(Meter))

# finally, reflectivities can be computed
let Au = Gold.init()
# note: currently densities of elements are not read by default. 
discard Au.plotReflectivity(19.32.g•cm⁻³,
                            θ = 0.5.°) # incidence angle of 0.5° from surface
# the procedure (currently, but that will be changed) returns the DF containing the reflectivity
#+end_src

See the API (docs will be published soon) / source file on how to compute individual properties
directly instead of only plotting the data.

** Data

The data used for the calculations are a mix of the dataset by NIST:

https://www.nist.gov/pml/x-ray-mass-attenuation-coefficients

and the Center for X-ray Optics ("Henke"):

https://henke.lbl.gov/optical_constants/

and specifically:

https://henke.lbl.gov/optical_constants/asf.html



** Resources

This section lists a few resources that can be useful for anyone who
wants to understand where the equations used here come from.

*** General
The X-ray data booklet:

https://xdb.lbl.gov/xdb-new.pdf

gives a good overview of all the topics required by this library (and
more of course!)

Note though that the book is short on maths and some the equations 
are hard to implement correctly (grazing angle reflectivity for
example).

And of course the main website of it:

https://xdb.lbl.gov/


In particular this PDF of "Useful formulas" is pretty handy:

https://xdb.lbl.gov/Section5/Sec_5-5.pdf

Further, NIST provides another database of scattering factors and
other numbers:
https://physics.nist.gov/PhysRefData/FFast/html/form.html

*** Other libraries for X-ray calculations

The DarpanX library provides similar functionality to this library,
with a focus on reflectivities for single and multi-layer mirrors (but
it's implemented in Python and rather slow):
- https://doi.org/10.1016/j.ascom.2020.100446
- https://github.com/biswajitmb/DarpanX

*** Reflectivity

The paper "Reflection of X-rays from a rough surface at extremely
small grazing angles"

https://doi.org/10.1364/OE.23.024220

was very useful in getting the basic reflectivity code working initially.

The following ~book contains an introduction about reflectivity,
including the derivation of the Fresnel equations (which provide the
basis to compute reflectivity and transmission), as well as

https://www.afc.asso.fr/images/reflecto2018/reflectie.pdf

Very short introduction to the topic from a lab course manual on X-ray
reflectometry at Uni Siegen:

https://www.hep.physik.uni-siegen.de/teaching/masterlab/manuals/XRR-2019-manual.pdf

Wikipedia on Fresnel equations:

https://en.wikipedia.org/wiki/Fresnel_equations

See also the paper about the DarpanX library, in particular the
appendix for an overview of the basic approach to compute
reflectivities.


*** Surface roughness

Mentions the origin of the dampening factor to Rayleigh & acoustic
waves.

https://www.classe.cornell.edu/~dms79/refl/XR-Roughness.html

also mentions Névot–Croce factors as a generalization of that.


Paper: "Influence of surface and interface roughness on X-ray and
extreme ultraviolet reflectance: A comparative numerical study"
- https://doi.org/10.1364/OSAC.422924
- https://opg.optica.org/osac/fulltext.cfm?uri=osac-4-5-1497&id=450674
seems to provide a good introduction.
