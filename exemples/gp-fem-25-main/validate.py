import numpy as np

ref_file = "data/UV_ref.txt"
sol_file = "data/UV.txt"

reference = np.loadtxt(ref_file, skiprows=1, delimiter=",")
solution = np.loadtxt(sol_file, skiprows=1, delimiter=",")

rel_err = np.linalg.norm(solution - reference) / np.linalg.norm(reference)
if rel_err > 1e-4:
    raise ValueError(f"Error is too large : {rel_err:.3e}")
else:
    print(f"Good job, relative error = {rel_err:.3e} < 1e-4")
