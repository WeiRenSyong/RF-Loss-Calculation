import numpy as np

# Constants
mu_0 = 4 * np.pi * 1e-7  # Vacuum permeability [H/m]

def RF_sheet_resistance(f_GHz, rho_ohm_over_m, mu_r):
    """
    Calculate the RF sheet resistance (R_RFSH) of a conductor.

    Input:
        f_Hz: Frequency in Hz.
        rho_ohm_over_m: Resistivity (ρ) of the conductor in ohm/m.
        mu_r: Relative permeability of the conductor.

    Output:
        R_RFSH: RF sheet resistance in ohms/square.
    """
    f_Hz = f_GHz * 1e9
    mu = mu_0 * mu_r  # Absolute permeability
    sigma = 1 / rho_ohm_over_m  # Conductivity
    R_RFSH = np.sqrt(np.pi * f_Hz * mu / sigma)
    return R_RFSH

# Example parameters
f_GHz = 1 
rho_copper = 1.67e-8  # Resistivity of copper in ohm-m
rho_aluminum = 2.65e-8  # Resistivity of aluminum in ohm-m
mu_r_copper = 0.999994
mu_r_aluminum = 1.00000065

# Compute sheet resistances
R_RFSH_copper = RF_sheet_resistance(f_GHz, rho_copper, mu_r_copper)
R_RFSH_aluminum = RF_sheet_resistance(f_GHz, rho_aluminum, mu_r_aluminum)

print(f"Copper R_RFSH at {f_GHz} GHz: {R_RFSH_copper:.4e} ohm/sq")
print(f"Aluminum R_RFSH at {f_GHz} GHz: {R_RFSH_aluminum:.4e} ohm/sq")

def coax_resistance_per_length(f_GHz, 
                                rho_inner, d_inner, mu_r_inner,
                                rho_outer, D_outer, mu_r_outer):
    """
    Estimate total resistance per unit length (R') of a coaxial cable due to metal conductivity.
    
    Input:
        f_Hz: Frequency in Hz.
        rho_inner: Resistivity of inner conductor [ohm·m].
        d_inner (float): Diameter of inner conductor [m].
        mu_r_inner (float): Relative permeability of inner conductor.
        rho_outer (float): Resistivity of outer conductor [ohm·m].
        D_outer (float): Inner diameter of outer conductor (i.e., jacket) [m].
        mu_r_outer (float): Relative permeability of outer conductor.
        
    Output:
        float: Total resistance per meter [ohm/m].
    """
    # Step 1: Get RF sheet resistances
    Rsh_inner = RF_sheet_resistance(f_GHz, rho_inner, mu_r_inner)
    Rsh_outer = RF_sheet_resistance(f_GHz, rho_outer, mu_r_outer)
    
    # Step 2: Compute resistance per unit length
    R_inner = Rsh_inner / (np.pi * d_inner)
    R_outer = Rsh_outer / (np.pi * D_outer)
    
    return R_inner, R_outer

# Example parameters
d_copper = 0.00102616
D_aluminum = 0.004572
R_inner, R_outer = coax_resistance_per_length(f_GHz, rho_copper, d_copper, mu_r_copper, rho_aluminum, D_aluminum, mu_r_aluminum)
print(f"Copper resistance per meter at {f_GHz} GHz: {R_inner:.4e} ohm/m")
print(f"Aluminum resistance per meter at {f_GHz} GHz: {R_outer:.4e} ohm/m")
