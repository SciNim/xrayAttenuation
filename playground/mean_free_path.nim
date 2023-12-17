import xrayAttenuation, datamancer


let niC = compound((N, 2))
let oC = compound((O, 2))
let arC = compound((Ar, 1)) # need Argon gas as a Compound
let co2C = compound((C, 1), (O, 2))
# define the gas mixture
let gm = initGasMixture(293.K, 1013.mbar, [(niC, 0.78084), (oC, 0.20946), (arC, 0.00934), (co2C, 0.000417)])

echo gm.absorptionLength(400.keV)

let lead = Lead.init()
echo lead.absorptionLength(11.34.g•cm⁻³, 600.keV)
