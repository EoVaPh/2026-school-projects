from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from core_phi_psi_functions import (read_alignment, read_angles, collect_core_angles, plot_scatter)


alignment_file = Path("MAFFT_families_atomseq_multiple_alignment\\family_004_aligned.txt")
angles_file = Path("all_phi_psi.txt")

alignment = read_alignment(alignment_file)
angles = read_angles(alignment, angles_file)
core, core_phis, core_psis = collect_core_angles(alignment, angles)


core_cos_phis = [[np.cos(phi) for phi in phis if phi is not None] for phis in core_phis]
core_cos_psis = [[np.cos(psi) for psi in psis if psi is not None] for psis in core_psis]
core_sin_phis = [[np.sin(phi) for phi in phis if phi is not None] for phis in core_phis]
core_sin_psis = [[np.sin(psi) for psi in psis if psi is not None] for psis in core_psis]

positions = np.arange(len(core))

plot_scatter(
    core_cos_phis,
    positions,
    "cos(φ)",
    "cos(φ) angles in the core",
    "core_phi_cos.png"
)
plot_scatter(
    core_cos_psis,
    positions,
    "cos(ψ)",
    "cos(ψ) angles in the core",
    "core_psi_cos.png"
)
plot_scatter(
    core_sin_phis,
    positions,
    "sin(φ)",
    "sin(φ) angles in the core",
    "core_phi_sin.png"
)
plot_scatter(
    core_sin_psis,
    positions,
    "sin(ψ)",
    "sin(ψ) angles in the core",
    "core_psi_sin.png"
)