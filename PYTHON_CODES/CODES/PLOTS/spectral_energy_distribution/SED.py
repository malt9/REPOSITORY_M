# SED.py
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.modeling.physical_models import BlackBody
from astropy.wcs import WCS
import astropy.constants as const
import astropy.units as u
from scipy.optimize import least_squares
from scipy.ndimage import map_coordinates


# =============================================================================
PROJECT_ROOT = Path(__file__).resolve().parent
DATA_DIR = PROJECT_ROOT / "data\\binned"
OUT_CSV = PROJECT_ROOT / "n113_sed_points_binned.csv"

# Única posición de interés: N113 (paper)
N113 = SkyCoord("05h13m17.40s", "-69d22m22.0s", frame="icrs")  # posicion de N113
N159W = SkyCoord("05h39m36.5s", "-69d45m35.0s", frame="icrs") 

# Extracción del valor de mapa (sin aperturas): bilinear es estable en bordes
EXTRACT_METHOD = "bilinear"  # "nearest" o "bilinear"

# Error fraccional por banda (paper no lo fija para la SED; deja esto explícito)
FRAC_ERR = 0.10


# =============================================================================
# Parámetros del paper (Apéndice A)
# =============================================================================
R_GDR = 300.0 * u.dimensionless_unscaled                       # LMC gas-to-dust ratio [file:1]
KAPPA_230 = 0.8 * (u.cm**2 / u.g)                              # κ_230GHz [file:1]
NU0 = 230.0 * u.GHz                                            # ν0 [file:1]
OMEGA_BEAM = 4.26e-8 * u.sr                                    # Ω para HPBW=40" [file:1]
BETA_FIXED = 1.96                                              # β fijo [file:1]
N_H_ATOMIC = 2.6e21 / (u.cm**2)                                # NH atómico [file:1]
M_H = const.m_p.to(u.g)                                        # m_H ~ masa del protón


@dataclass(frozen=True)
class BandSpec:
    name: str
    wav: u.Quantity
    pattern: re.Pattern


BANDS = [
    BandSpec("PACS100", 100.0 * u.um, re.compile(r"lmc_both_pacs_100_.*\.fits$", re.IGNORECASE)),
    BandSpec("PACS160", 160.0 * u.um, re.compile(r"lmc_both_pacs_160_.*\.fits$", re.IGNORECASE)),
    BandSpec("SPIRE250", 250.0 * u.um, re.compile(r"lmc_both_spire_250_.*\.fits$", re.IGNORECASE)),
    BandSpec("SPIRE350", 350.0 * u.um, re.compile(r"lmc_both_spire_350_.*\.fits$", re.IGNORECASE)),
    BandSpec("SPIRE500", 500.0 * u.um, re.compile(r"lmc_both_spire_500_.*\.fits$", re.IGNORECASE)),
]


# =============================================================================
# Modelo físico (Apéndice A)
# =============================================================================
def kappa_nu(nu: u.Quantity, beta: float) -> u.Quantity:
    """κ_ν = κ_0 (ν/ν0)^β   (cm^2/g)."""
    nu = nu.to(u.Hz)
    return KAPPA_230 * (nu / NU0.to(u.Hz)) ** beta


def tau_nu(nu: u.Quantity, Ntot: u.Quantity, beta: float) -> u.Quantity:
    """τ_ν = Ntot κ_ν m_H / rGDR."""
    Ntot = Ntot.to(1 / u.cm**2)
    return (Ntot * kappa_nu(nu, beta) * M_H / R_GDR).to(u.dimensionless_unscaled)


def model_snu(nu: u.Quantity, T: u.Quantity, Ntot: u.Quantity, beta: float) -> u.Quantity:
    """S_ν = Ω B_ν(T) (1 - exp(-τ_ν))   en Jy."""
    bb = BlackBody(temperature=T)
    Bnu = bb(nu).to(u.Jy / u.sr)
    tau = tau_nu(nu, Ntot, beta)
    return (OMEGA_BEAM * Bnu * (1.0 - np.exp(-tau.value))).to(u.Jy)


# =============================================================================
# Datos: 
#   - localizar FITS 
#   - extraer valores
# =============================================================================
def find_file_for_band(spec: BandSpec) -> Path:
    candidates = [p for p in DATA_DIR.iterdir() if p.is_file() and spec.pattern.search(p.name)]
    if not candidates:
        raise FileNotFoundError(f"No se encontró FITS para {spec.name} en {DATA_DIR}")
    if len(candidates) > 1:
        msg = "\n".join(f" - {p.name}" for p in candidates)
        raise RuntimeError(f"Más de un FITS coincide para {spec.name}:\n{msg}")
    return candidates[0]


def read_map_value(fits_path: Path, coord: SkyCoord, method: str) -> u.Quantity:
    """Devuelve I_ν en MJy/sr en la coordenada dada."""
    with fits.open(fits_path) as hdul:
        hdu = hdul[0]
        data = np.asarray(hdu.data, dtype=float)
        wcs = WCS(hdu.header)

    if data.ndim != 2:
        raise ValueError(f"{fits_path.name}: se esperaba imagen 2D, data.ndim={data.ndim}")

    xpix, ypix = wcs.world_to_pixel_values(coord.ra.deg, coord.dec.deg)

    if method == "nearest":
        xi, yi = int(np.round(xpix)), int(np.round(ypix))
        val = data[yi, xi]
    elif method == "bilinear":
        val = map_coordinates(data, [[ypix], [xpix]], order=1, mode="nearest")[0]
    else:
        raise ValueError("method debe ser 'nearest' o 'bilinear'")

    return val * (u.MJy / u.sr)


def intensity_to_flux(I: u.Quantity) -> u.Quantity:
    """S_ν = I_ν Ω."""
    return (I.to(u.Jy / u.sr) * OMEGA_BEAM).to(u.Jy)


# =============================================================================
# Ajuste
# =============================================================================
def fit_sed(wavs: u.Quantity, S: u.Quantity, Serr: u.Quantity, beta_mode: str):
    nu = (const.c / wavs.to(u.m)).to(u.Hz)
    y = S.to_value(u.Jy)
    yerr = Serr.to_value(u.Jy)

    if beta_mode == "fixed":
        beta = float(BETA_FIXED)

        # params: [T_K, log10(Ntot/cm^-2)]
        p0 = np.array([20.0, 22.8])
        lo = np.array([3.0, 20.0])
        hi = np.array([150.0, 25.0])

        def resid(p):
            T = p[0] * u.K
            N = (10.0 ** p[1]) * (1 / u.cm**2)
            m = model_snu(nu, T, N, beta).to_value(u.Jy)
            return (y - m) / yerr

        res = least_squares(resid, p0, bounds=(lo, hi))
        return {"T": res.x[0] * u.K, "Ntot": (10.0 ** res.x[1]) / (u.cm**2), "beta": beta,
                "success": res.success, "cost": res.cost}

    if beta_mode == "free":
        # params: [T_K, log10(Ntot/cm^-2), beta]
        p0 = np.array([30.0, 22.7, 1.2])
        lo = np.array([3.0, 20.0, -1.0])
        hi = np.array([150.0, 25.0, 4.0])

        def resid(p):
            T = p[0] * u.K
            N = (10.0 ** p[1]) * (1 / u.cm**2)
            beta = float(p[2])
            m = model_snu(nu, T, N, beta).to_value(u.Jy)
            return (y - m) / yerr

        res = least_squares(resid, p0, bounds=(lo, hi))
        return {"T": res.x[0] * u.K, "Ntot": (10.0 ** res.x[1]) / (u.cm**2), "beta": float(res.x[2]),
                "success": res.success, "cost": res.cost}

    raise ValueError("beta_mode debe ser 'fixed' o 'free'")


# =============================================================================
# Main
# =============================================================================
def main():
    rows = []
    for spec in BANDS:
        path = find_file_for_band(spec)
        I = read_map_value(path, N113, method=EXTRACT_METHOD)      # MJy/sr
        S = intensity_to_flux(I)                                   # Jy

        rows.append({
            "band": spec.name,
            "wav_um": spec.wav.to_value(u.um),
            "file": path.name,
            "I_MJy_sr": I.to_value(u.MJy / u.sr),
            "S_Jy": S.to_value(u.Jy),
        })

    df = pd.DataFrame(rows).sort_values("wav_um")
    df["Serr_Jy"] = FRAC_ERR * df["S_Jy"]
    df.to_csv(OUT_CSV, index=False)

    wav = df["wav_um"].to_numpy() * u.um
    S = df["S_Jy"].to_numpy() * u.Jy
    Serr = df["Serr_Jy"].to_numpy() * u.Jy

    fit_fixed = fit_sed(wav, S, Serr, beta_mode="fixed")
    fit_free = fit_sed(wav, S, Serr, beta_mode="free")

    def report(tag: str, fit: dict):
        Ntot = fit["Ntot"].to(1 / u.cm**2)
        NH2 = ((Ntot - N_H_ATOMIC) / 2.0).to(1 / u.cm**2)
        if NH2.value < 0:
            NH2 = 0 * (1 / u.cm**2)

        print(f"\n[{tag}]")
        print(f"T_d [K]        = {fit['T']:.3f}")
        print(f"beta           = {fit['beta']:.3f}")
        print(f"Ntot(H) [cm^-2]= {Ntot:.3e}")
        print(f"N(H2)  [cm^-2] = {NH2:.3e}")
        print(f"CSV            = {OUT_CSV.name}")
        print(f"success        = {fit['success']}, cost={fit['cost']:.3g}")

    report("beta fixed (paper: beta=1.96)", fit_fixed)             # [file:1]
    report("beta free  (paper: better fit)", fit_free)             # [file:1]


if __name__ == "__main__":
    main()
