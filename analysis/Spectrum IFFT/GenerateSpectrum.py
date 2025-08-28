# Generate a two-column spectrum text file:
# Column 1: wavelength in nm (700 to 900 nm, step 1 nm)
# Column 2: intensity (Gaussian centered at 815 nm, FWHM = 35.8 nm) with a small amount of noise

import numpy as np

# Parameters
center_nm = 815.0
fwhm_nm = 35.8
wavelengths = np.arange(700.0, 901.0, 1.0)  # 700..900 nm inclusive, 1 nm step

# Gaussian intensity profile: I = exp(-4*ln(2) * ((λ - λ0)/FWHM)^2)
gaussian = np.exp(-4*np.log(2) * ((wavelengths - center_nm) / fwhm_nm)**2)

# Add small nonnegative noise (simulate measurement baseline and noise floor)
rng = np.random.default_rng(42)
intensity = np.clip(gaussian, 0, None)

# Save as plain text with two columns, space-separated, no header (compatible with MATLAB load())
out_path = "gaussian spectrum.txt"
np.savetxt(out_path, np.column_stack([wavelengths, intensity]), fmt="%.6f")

out_path
