m_dot_fuel =0.0393
v_fuel= 29.5
rho_fuel= 832.59
rho_ox= 1450
A_fuel= m_dot_fuel/(rho_fuel*v_fuel)
print(A_fuel)

#TMR =( m_dot_fuel * v_fuel )/(m_dot_ox * v_ox)

m_dot_ox=5*m_dot_fuel
TMR_values = [0.75,1,1.25,1.5] # Renamed to TMR_values for clarity
v_ox_list = [m_dot_fuel/m_dot_ox * v_fuel / tmr for tmr in TMR_values] # Calculate v_ox for each TMR value
print(v_ox_list)
A_ox_list = [m_dot_ox/(rho_ox*v_ox_val) for v_ox_val in v_ox_list] # Calculate A_ox for each v_ox value, corrected rho to rho_ox
print(A_ox_list)


ID = 9e-3 #mm
# Calculate the fixed part of the area once
fixed_area_part = 3.14 * 0.25 * (ID**2)

# Calculate 'area' for each element in A_ox_list
area_list = [fixed_area_part + a_ox_val for a_ox_val in A_ox_list]

# Calculate 'OD' for each corresponding 'area'
OD_list = [(4 * current_area / 3.14)**0.5 * 1000 for current_area in area_list]##converting to mm

print("\n",OD_list)
#print(area_list)
