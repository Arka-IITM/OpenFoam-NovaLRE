# ============================================================
#  orifice_sizing.m — Injector Hole Area and Geometry
#  Logic: TMR and BR are inputs, dP is output
# ============================================================

# The previous orifice sizing was giving a really small orifice size so I
# tried using the full bernoulli instead of only the orifice velocity.
# Iterating to find A_fuel.


# Observation- No significant change in fuel area... Means what I calculated
# isnt incorrect..but I still doubt my calculations

# These variables will be overridden by the params cell if it's run before this cell
# mdot_fuel=0.0393;
# mdot_ox= 5 * mdot_fuel;
# dP_max=1000000

# ---- DISCHARGE COEFFICIENT ---------------------------------
# Cd_fuel  = 0.35;
# Cd_ox    = 0.35;

# Using variables defined in the params cell

# ---- STEP 1: FUEL TOTAL AREA FROM ORIFICE EQUATION --------
# A = mdot / (Cd * sqrt(2 * dP * rho))
# Using max available dP
A_fuel   = mdot_fuel / (Cd_fuel * math.sqrt(2 * dP_max * rho_fuel));

# ---- STEP 2: OX TOTAL AREA FROM TMR -----------------------
# TMR = (mdot_fuel^2 * rho_ox * A_ox) / (mdot_ox^2 * rho_fuel * A_fuel)
# Solving for A_ox:a
A_ox     = TMR_target * (mdot_ox**2 * rho_fuel * A_fuel) / (mdot_fuel**2 * rho_ox);

# ---- STEP 3: VELOCITIES ------------------------------------
V_fuel   = mdot_fuel / (rho_fuel * A_fuel);
V_ox     = mdot_ox   / (rho_ox   * A_ox);

# ---- STEP 4: VERIFY TMR ------------------------------------
TMR_actual = (mdot_fuel * V_fuel) / (mdot_ox * V_ox);

# ---- STEP 5: RESULTING dP ----------------------------------
dP_fuel  = 0.5 * rho_fuel * (V_fuel / Cd_fuel)**2;
dP_ox    = 0.5 * rho_ox   * (V_ox   / Cd_ox)**2;

# ---- PRINT RESULTS -----------------------------------------
print('\n--- ORIFICE SIZING RESULTS ---\n')
print('FUEL SIDE:')
print(f'  A_fuel total     : {A_fuel*1e6:.4f} mm2')
print(f'  V_fuel           : {V_fuel:.2f} m/s')
print(f'  dP_fuel          : {dP_fuel/1e5:.4f} bar')
print('\nOXIDISER SIDE:')
print(f'  A_ox total       : {A_ox*1e6:.4f} mm2')
print(f'  V_ox             : {V_ox:.2f} m/s')
print(f'  dP_ox            : {dP_ox/1e5:.4f} bar')
print('\nVERIFICATION:')
print(f'  TMR actual       : {TMR_actual:.4f} (target: {TMR_target:.2f})')

# ---- STEP 6: OX ANNULAR GAP GEOMETRY -----------------------
# Annular area: A_ox = pi/4 * (D_outer^2 - D_pintle^2)
# Solving for D_outer:
D_outer_annulus = math.sqrt((4 * A_ox / math.pi) + D_pintle**2);
gap_ox          = (D_outer_annulus - D_pintle) / 2;

# ---- CHECK GAP FITS WITHIN CHAMBER -------------------------
gap_available   = (D_chamber - D_pintle) / 2;
gap_remaining   = gap_available - gap_ox;

print('\nOX ANNULAR GAP:')
print(f'  D_outer annulus  : {D_outer_annulus*1e3:.4f} mm')
print(f'  Annular gap      : {gap_ox*1e3:.4f} mm')
print(f'  Available gap    : {gap_available*1e3:.4f} mm')
print(f'  Remaining margin : {gap_remaining*1e3:.4f} mm')

if gap_remaining < 0:
    print('  WARNING: Annular gap EXCEEDS chamber wall')
else:
    print('  Gap check        : OK')

# ---- STEP 8: FUEL SLOT SIZING ------------------------------
# Slot length from BR
l_slot       = (BR_slot * math.pi * D_pintle) / n_slots;

# Bridge width from circumference constraint
b_bridge     = (math.pi * D_pintle / n_slots) - l_slot;

# Slot height from flow area (single row)
# A_fuel = n_slots * w * l_slot
w_slot       = A_fuel / (n_slots * l_slot);

# Verify BR
BR_slot_actual = (n_slots * l_slot) / (math.pi * D_pintle);

print('\nFUEL SLOT SIZING:')
print(f'  Number of slots  : {n_slots:d}')
print(f'  Slot length      : {l_slot*1e3:.4f} mm')
print(f'  Slot height      : {w_slot*1e3:.4f} mm')
print(f'  Bridge width     : {b_bridge*1e3:.4f} mm')
print(f'  BR actual        : {BR_slot_actual:.4f}')
print(f'  BR target        : {BR_slot:.4f}')

