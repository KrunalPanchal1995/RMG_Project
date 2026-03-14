# input.py — RMG input file for NH3 oxidation (air) mechanism generation

database(
    thermoLibraries=['primaryThermoLibrary'],
    reactionLibraries=[],
    seedMechanisms=[],
    kineticsDepositories=['training'],
    kineticsFamilies='default',
    kineticsEstimator='rate rules',
)

# --- Species definitions ---
species(
    label='NH3',
    reactive=True,
    structure=SMILES("N"),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)

# Inert bath gas for air
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N"),
)

# --- Reactor conditions ---
# Choose conditions relevant to your problem (T/P/mixture).
# The mole fractions below are an example "lean" mixture.
simpleReactor(
    temperature=(1100, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        "NH3": 0.03,
        "O2": 0.18,
        "N2": 0.79,
    },
    terminationTime=(0.02, 's'),
)

# (Optional) Add more reactors to broaden the chemistry captured
# Example: slightly lower T, higher P (relevant for engines/gas turbines)
# simpleReactor(
#     temperature=(900, 'K'),
#     pressure=(10.0, 'bar'),
#     initialMoleFractions={
#         "NH3": 0.03,
#         "O2": 0.18,
#         "N2": 0.79,
#     },
#     terminationTime=(0.05, 's'),
# )

simulator(
    atol=1e-16,
    rtol=1e-8,
)

# --- Model construction settings ---
# For NH3, you may need smaller toleranceMoveToCore (e.g., 0.05 or 0.01)
# if you want more complete nitrogen chemistry (NHx/NOx paths).
model(
    toleranceMoveToCore=0.1,
    toleranceKeepInEdge=0.0,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=30000,
)

options(
    units='si',
    generateSeedEachIteration=False,
    generateOutputHTML=True,
    saveSimulationProfiles=True,
    verboseComments=False,
)
