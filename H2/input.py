# input.py — RMG input file for H2 oxidation mechanism generation

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
    label='H2',
    reactive=True,
    structure=SMILES("[H][H]"),
)

species(
    label='O2',
    reactive=True,
    structure=SMILES("[O][O]"),
)

# Inert bath gas (optional but recommended for realistic conditions)
species(
    label='N2',
    reactive=False,
    structure=SMILES("N#N"),
)

# --- Reactor conditions ---
# A constant T/P batch reactor. Choose conditions relevant to your use case.
simpleReactor(
    temperature=(1000, 'K'),
    pressure=(1.0, 'bar'),
    initialMoleFractions={
        "H2": 0.05,
        "O2": 0.21,
        "N2": 0.74,
    },
    terminationTime=(0.02, 's'),
)

# (Optional) Add additional reactors to broaden the generated mechanism
# Example: lower temperature condition to encourage inclusion of low-T pathways
# simpleReactor(
#     temperature=(800, 'K'),
#     pressure=(1.0, 'bar'),
#     initialMoleFractions={
#         "H2": 0.05,
#         "O2": 0.21,
#         "N2": 0.74,
#     },
#     terminationTime=(0.05, 's'),
# )

simulator(
    atol=1e-16,
    rtol=1e-8,
)

# --- Model construction settings ---
# These tolerances control mechanism growth:
# smaller toleranceMoveToCore => larger mechanism, more complete chemistry
model(
    toleranceMoveToCore=0.1,
    toleranceKeepInEdge=0.0,
    toleranceInterruptSimulation=1.0,
    maximumEdgeSpecies=20000,
)

options(
    units='si',
    generateSeedEachIteration=False,
    generateOutputHTML=True,
    saveSimulationProfiles=True,
    verboseComments=False,
)
