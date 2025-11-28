h1 = 500e3
h2 = 35789e3
deltaV1, deltaV2 = earth.get_Hohmann_transfer_deltaV(h1, h2)
mass2 = simlib.Engine(None, simlib.ariane_stage2).compute_propellant_mass(deltaV2, dry_mass=1e3+4.54e3)
mass1 = simlib.Engine(None, simlib.ariane_stage2).compute_propellant_mass(deltaV1, dry_mass=1e3+4.54e3 + mass2)

# LIFTOFF
STAGE_1_PROPELLANT_MASS = 152e3
ROTATION_ANGLE = -75

# HOHMANN TRANSFER
STAGE_2_TOTAL_MASS = 12e3
STAGE_2_FIRST_BURN_MASS = 7324.5
STAGE_2_SECOND_BURN_MASS = 1600