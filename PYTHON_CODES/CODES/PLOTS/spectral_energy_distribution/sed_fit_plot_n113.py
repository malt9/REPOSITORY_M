
# Uso:
#   python SED.py n113_sed_points_un / binned.csv


import sys
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

# -----------------------------
# Constantes / parámetros (Gong+2024, Apéndice A)
# -----------------------------
h = 6.62607015e-27       # erg s
kB = 1.380649e-16        # erg K^-1
c = 2.99792458e10        # cm s^-1
mH = 1.6735575e-24       # g

Omega = 4.26e-8          # sr (HPBW 40")
rGDR = 300.0             # gas-to-dust ratio (LMC)
kappa0 = 0.8             # cm^2 g^-1 at nu0
nu0 = 230e9              # Hz

def Bnu_cgs(nu_hz, T):
    """Planck B_nu en cgs: erg s^-1 cm^-2 Hz^-1 sr^-1"""
    x = h * nu_hz / (kB * T)
    return (2.0 * h * nu_hz**3 / c**2) / np.expm1(x)

def Snu_model_Jy(wav_um, T, NH, beta):
    """Devuelve S_nu en Jy (beam 40") para wav_um (micrones)."""
    wav_cm = wav_um * 1e-4
    nu_hz = c / wav_cm
    kappa = kappa0 * (nu_hz / nu0)**beta
    tau = NH * mH * kappa / rGDR
    I_nu = Bnu_cgs(nu_hz, T) * (1.0 - np.exp(-tau))  # cgs
    I_nu_Jy_sr = I_nu / 1e-23                        # Jy sr^-1
    return Omega * I_nu_Jy_sr                        # Jy

def model_beta_free(wav_um, T, NH, beta):
    return Snu_model_Jy(wav_um, T, NH, beta)

def model_beta_fixed(beta):
    def _m(wav_um, T, NH):
        return Snu_model_Jy(wav_um, T, NH, beta)
    return _m

def fit_fixed_beta(wav_um, S_Jy, Serr_Jy, beta, p0=(30.0, 6e22)):
    popt, pcov = curve_fit(
        model_beta_fixed(beta),
        wav_um, S_Jy,
        sigma=Serr_Jy,
        p0=p0,
        bounds=([5.0, 1e20], [150.0, 5e24]),
        absolute_sigma=True,
        maxfev=200000
    )
    return popt, pcov  # (T, NH)

def fit_free_beta(wav_um, S_Jy, Serr_Jy, p0=(35.0, 6e22, 1.2)):
    popt, pcov = curve_fit(
        model_beta_free,
        wav_um, S_Jy,
        sigma=Serr_Jy,
        p0=p0,
        bounds=([5.0, 1e20, -1.0], [150.0, 5e24, 4.0]),
        absolute_sigma=True,
        maxfev=300000
    )
    return popt, pcov  # (T, NH, beta)

def fmt_NH(NH):
    exp = int(np.floor(np.log10(NH)))
    mant = NH / 10**exp
    return f"{mant:.2g}×10^{exp:d}"

def main(csv_path):
    df = pd.read_csv(csv_path)
    wav_um = df["wav_um"].to_numpy(float)
    S = df["S_Jy"].to_numpy(float)
    Serr = df["Serr_Jy"].to_numpy(float)

    # Tests 
    assert np.all(np.isfinite(wav_um)) and np.all(wav_um > 0)
    assert np.all(np.isfinite(S)) and np.all(S > 0)
    assert np.all(np.isfinite(Serr)) and np.all(Serr > 0)

    # Ajuste: beta fijo = 1.96 (usado en Gong et al. 2024)
    beta_fix = 1.96
    (T_fix, NH_fix), _ = fit_fixed_beta(wav_um, S, Serr, beta_fix, p0=(25.0, 6e22))

    # Ajuste: beta libre
    (T_free, NH_free, beta_free), _ = fit_free_beta(wav_um, S, Serr, p0=(35.0, 6e22, 1.2))

    # Diagnóstico: beta=1.18±0.15 (obtenido de Gong. et al.2024)
    beta0 = 1.18
    dbeta = 0.15
    beta_lo = beta0 - dbeta
    beta_hi = beta0 + dbeta

    (T_b0, NH_b0), _ = fit_fixed_beta(wav_um, S, Serr, beta0, p0=(35.0, 6e22))
    (T_blo, NH_blo), _ = fit_fixed_beta(wav_um, S, Serr, beta_lo, p0=(38.0, 6e22))
    (T_bhi, NH_bhi), _ = fit_fixed_beta(wav_um, S, Serr, beta_hi, p0=(33.0, 6e22))

    # Grid para curvas
    wgrid = np.logspace(np.log10(10.0), np.log10(2300.0), 600)

    y_fix = Snu_model_Jy(wgrid, T_fix, NH_fix, beta_fix)
    y_free = Snu_model_Jy(wgrid, T_free, NH_free, beta_free)

    # Banda del diagnóstico: envolvente entre beta_lo y beta_hi
    y_lo = Snu_model_Jy(wgrid, T_blo, NH_blo, beta_lo)
    y_hi = Snu_model_Jy(wgrid, T_bhi, NH_bhi, beta_hi)
    y_band_min = np.minimum(y_lo, y_hi)
    y_band_max = np.maximum(y_lo, y_hi)

    # Curva central (beta0)
    y_b0 = Snu_model_Jy(wgrid, T_b0, NH_b0, beta0)

    # ---- Plot  ----
    plt.figure(figsize=(6.2, 4.8), dpi=160)

    plt.errorbar(
        wav_um, S, yerr=Serr,
        fmt="o", color="black", ecolor="black",
        elinewidth=1.0, capsize=2.5, markersize=4.5
    )

    plt.plot(wgrid, y_free, color="#1f77b4", lw=2.0,
             label=f"$T_d$={T_free:.0f} K  $N$(H)={fmt_NH(NH_free)} cm$^{{-2}}$  $\\beta$={beta_free:.2f}")
    plt.plot(wgrid, y_fix, color="#ff7f0e", lw=2.0,
             label=f"$T_d$={T_fix:.0f} K  $N$(H)={fmt_NH(NH_fix)} cm$^{{-2}}$  $\\beta$={beta_fix:.2f}")

    plt.fill_between(wgrid, y_band_min, y_band_max, color="#2ca02c", alpha=0.25,
                     label=f"diag: $\\beta$={beta0:.2f}±{dbeta:.2f}")
    plt.plot(wgrid, y_b0, color="#2ca02c", lw=2.0,
             label=f"diag: $T_d$={T_b0:.0f} K  $N$(H)={fmt_NH(NH_b0)} cm$^{{-2}}$")

    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Wavelength (µm)")
    plt.ylabel("$S_\\nu$ (Jy)")
    plt.xlim(30, 2100)
    plt.ylim(1e-1, 10e2)

    plt.legend(frameon=False, fontsize=8, loc="upper right")
    plt.tight_layout()
    plt.savefig("n113_SED_fit_unbinned_pixel.pdf", dpi=300)
    plt.show()

    print("=== Fit summary ===")
    print(f"beta fixed {beta_fix:.2f}: Td={T_fix:.3f} K, NH={NH_fix:.6e} cm^-2")
    print(f"beta free: Td={T_free:.3f} K, NH={NH_free:.6e} cm^-2, beta={beta_free:.4f}")
    print(f"diag beta {beta0:.2f}: Td={T_b0:.3f} K, NH={NH_b0:.6e} cm^-2")
    print(f"diag beta {beta_lo:.2f}: Td={T_blo:.3f} K, NH={NH_blo:.6e} cm^-2")
    print(f"diag beta {beta_hi:.2f}: Td={T_bhi:.3f} K, NH={NH_bhi:.6e} cm^-2")
    print("Saved: n113_SED_fit_unbinned_pixel.pdf")

if __name__ == "__main__":
    csv_path = sys.argv[1] if len(sys.argv) > 1 else "n113_sed_points_unbinned.csv"
    main(csv_path)
