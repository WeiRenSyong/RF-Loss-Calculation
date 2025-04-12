import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
mu_0 = 4 * np.pi * 1e-7  # Vacuum permeability [H/m]


# --- Functions ---
def RF_sheet_resistance(f_GHz, rho_ohm_m, mu_r):
    """
    Calculate RF sheet resistance (R_RFSH) of a conductor.
    Args:
        f_GHz (float): Frequency in GHz.
        rho_ohm_m (float): Resistivity in ohm·m.
        mu_r (float): Relative permeability.
    Returns:
        float: Sheet resistance in ohm/sq.
    """
    f_Hz = f_GHz * 1e9
    mu = mu_0 * mu_r
    sigma = 1 / rho_ohm_m
    return np.sqrt(np.pi * f_Hz * mu / sigma)


def coax_resistance_per_length(f_GHz, rho_inner, d_inner, mu_r_inner,
                                rho_outer, D_outer, mu_r_outer):
    """
    Total resistance per unit length of a coaxial cable.
    Returns:
        R_inner, R_outer (float): Inner and outer conductor resistance in ohm/m.
    """
    Rsh_inner = RF_sheet_resistance(f_GHz, rho_inner, mu_r_inner)
    Rsh_outer = RF_sheet_resistance(f_GHz, rho_outer, mu_r_outer)
    R_inner = Rsh_inner / (np.pi * d_inner)
    R_outer = Rsh_outer / (np.pi * D_outer)
    return R_inner, R_outer


def coax_characteristic_impedance(d_inner, D_outer, epsilon_r):
    """
    Characteristic impedance of a coaxial cable.
    """
    return 138 / np.sqrt(epsilon_r) * np.log(D_outer / d_inner)


def coax_loss_per_length(f_GHz, rho_inner, d_inner, mu_r_inner,
                         rho_outer, D_outer, mu_r_outer,
                         epsilon_r):
    """
    Compute loss per unit length in dB/m.
    """
    R_inner, R_outer = coax_resistance_per_length(f_GHz, rho_inner, d_inner, mu_r_inner,
                                                  rho_outer, D_outer, mu_r_outer)
    Z0 = coax_characteristic_impedance(d_inner, D_outer, epsilon_r)
    loss_dB_per_m = 8.686 * (R_inner + R_outer) / (2 * Z0)
    return loss_dB_per_m


def plot_loss_vs_frequency(freq_range_GHz, rho_inner, d_inner, mu_r_inner,
                           rho_outer, D_outer, mu_r_outer,
                           epsilon_r):
    """
    Plot loss per length (dB/m) vs frequency.
    """
    freqs = np.linspace(*freq_range_GHz, 500)
    losses = [coax_loss_per_length(f, rho_inner, d_inner, mu_r_inner,
                                   rho_outer, D_outer, mu_r_outer,
                                   epsilon_r) for f in freqs]

    plt.figure(figsize=(8, 5))
    plt.plot(freqs, losses, label='Loss per meter', color='darkblue')
    plt.xlabel("Frequency (GHz)")
    plt.ylabel("Loss (dB/m)")
    plt.title("Coaxial Cable Loss vs Frequency")
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()


# --- Example Inputs ---
# Material properties
rho_copper = 1.67e-8       # [ohm·m]
rho_aluminum = 2.65e-8     # [ohm·m]
mu_r_copper = 0.999994
mu_r_aluminum = 1.00000065

# Geometry
d_inner = 0.00102616       # [m]
D_outer = 0.004572         # [m]
epsilon_r = 1

# --- Example: Plot loss vs frequency ---
freq_range = (0.1, 3)  # GHz, from 10 MHz to 10 GHz
plot_loss_vs_frequency(freq_range,
                       rho_copper, d_inner, mu_r_copper,
                       rho_aluminum, D_outer, mu_r_aluminum,
                       epsilon_r)