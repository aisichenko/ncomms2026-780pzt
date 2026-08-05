"""Resonator Q-fitting utilities for the NCOMMS data repository."""
# Credit: Portions of this code are adapted from the pyLaserNoise repository
# and from Kaikai Liu's pyphotonics simulations (pyphotonics_sims).


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.constants import speed_of_light as c
from scipy.optimize import curve_fit
from scipy.signal import find_peaks



def linear_fit(x, y):
    """Least-squares linear fit y = m*x + b. Returns (m, b)."""
    m, b = np.polyfit(x, y, 1)
    return m, b


def TAddThru(dw, r_in, r_ex, w0=0, level=1):
    """Add-thru resonator cavity transmission.

    Args:
        dw: frequency detuning in [MHz]
        r_in: intrinsic loss in [MHz]
        r_ex: external loss in [MHz]
        w0: resonance offset in [MHz]
        level: transmission scale

    Returns:
        T: cavity transmission
    """
    dw = dw - w0
    F = (1j * dw + (r_in - r_ex) / 2) / (1j * dw + (r_in + r_ex) / 2)
    return level * abs(F) ** 2


def TAddThruSplit(dw, r_in, r_ex, g, w0=0, level=1):
    """Add-thru resonator with CW-CCW coupling (split resonance).

    Args:
        dw: detuning in [MHz]
        r_in: intrinsic loss in [MHz]
        r_ex: external loss in [MHz]
        g: CW-CCW coupling rate in [MHz]
        w0: resonance offset in [MHz]
        level: transmission scale

    Returns:
        T: cavity transmission
    """
    dw = dw - w0
    F = 1 - r_ex * (1j * dw + (r_in + r_ex) / 2) / ((1j * dw + (r_in + r_ex) / 2) ** 2 + g ** 2)
    return level * abs(F) ** 2


def moving_average(x, w):
    return np.convolve(x, np.ones(w), "valid") / w


def crop_rescale_res_mzi(path, trans_col=1, mzi_col=2, ramp_col=3, X_window=0, Nmean=1, mzi_debug=True):
    """Load resonance/MZI/ramp CSV and optionally crop around the resonance.

    Args:
        path: path to the csv file
        X_window: crop to X times the resonance linewidth (0 = no crop)
        Nmean: points for moving average
        mzi_debug: plot loaded traces

    Returns:
        y_res: normalized resonance data
        y_mzi: normalized MZI data
    """
    df = pd.read_csv(path)
    df.dropna(inplace=True)

    required_columns = ["trans_pd", "ramp", "mzi_pd"]
    df_check = set(required_columns).issubset(list(df.columns))

    if not df_check:
        df = pd.read_csv(path, skiprows=10)
        df.dropna(inplace=True)
        cols = df.columns
        df["trans_pd"] = df[cols[trans_col]]
        df["mzi_pd"] = df[cols[mzi_col]]
        df["ramp"] = df[cols[ramp_col]]
        df = df[["trans_pd", "mzi_pd", "ramp"]]

    y_res = df["trans_pd"].to_numpy()
    y_res = moving_average(y_res, Nmean)
    y_res = y_res / np.max(y_res)

    y_mzi = df["mzi_pd"].to_numpy()
    y_mzi = moving_average(y_mzi, Nmean)
    y_mzi = y_mzi - np.min(y_mzi)
    y_mzi = y_mzi / np.max(y_mzi) * 0.2

    y_ramp = df["ramp"].to_numpy()
    y_ramp = moving_average(y_ramp, Nmean)
    y_ramp = y_ramp / np.max(y_ramp) * 0.5

    if X_window != 0:
        idx_res = np.argmin(y_res)
        y_fwhm = 0.5 * (np.max(y_res) + np.min(y_res))
        idx_fwhm = np.argmin(np.abs(y_res[0:idx_res] - y_fwhm))
        idx_start = int(idx_res - int(2 * X_window * (idx_res - idx_fwhm)))
        idx_end = int(idx_res + int(2 * X_window * (idx_res - idx_fwhm)))

        if idx_start < 0:
            idx_start = 1
        if idx_end > len(y_res) - 1:
            idx_end = len(y_res) - 2

        y_res = y_res[idx_start:idx_end]
        y_mzi = y_mzi[idx_start:idx_end]
        y_ramp = y_ramp[idx_start:idx_end]

    if mzi_debug:
        plt.figure()
        plt.plot(y_res, label="Resonance")
        plt.plot(y_mzi, label="MZI")
        plt.plot(y_ramp, label="Ramp")
        plt.xlabel("X")
        plt.ylabel("Transmission (a.u.)")
        plt.title("Crop and rescale data")
        plt.legend()
        plt.tight_layout()
        plt.show()
    return y_res, y_mzi


def calibrate_mzi(y_mzi, fsr_mzi=5.87, Nmean=1, Nbetween=0, peaks_debug=True):
    """Calibrate frequency axis from MZI fringe peaks.

    Args:
        y_mzi: mzi data
        fsr_mzi: FSR of the MZI fringes in MHz
        Nmean: points for moving average
        Nbetween: minimum distance between peaks
        peaks_debug: plot peak finding

    Returns:
        f: frequency detuning [MHz]
    """
    y_mzi_1 = moving_average(y_mzi, Nmean)

    if Nbetween == 0:
        ind_p, _ = find_peaks(y_mzi_1, height=0.1)
    else:
        ind_p, _ = find_peaks(y_mzi_1, distance=Nbetween, height=0.1)

    f = np.linspace(1, len(y_mzi), len(y_mzi))
    f1 = np.linspace(1, len(y_mzi_1), len(y_mzi_1))
    f = fsr_mzi / np.mean(np.diff(ind_p)) * f
    f1 = fsr_mzi / np.mean(np.diff(ind_p)) * f1

    if peaks_debug:
        plt.figure()
        plt.subplot(2, 1, 1)
        plt.plot(y_mzi_1)
        plt.plot(ind_p, y_mzi_1[ind_p], "x")
        plt.xlabel("Index")
        plt.subplot(2, 1, 2)
        if len(ind_p) < 50:
            plt.plot(f1, y_mzi_1)
            plt.plot(f1[ind_p], y_mzi_1[ind_p], "x")
        else:
            n_ind = len(ind_p)
            ind_center = int(n_ind / 2)
            plt.plot(f1, y_mzi_1)
            plt.plot(
                f1[ind_p[ind_center - 10 : ind_center + 10]],
                y_mzi_1[ind_p[ind_center - 10 : ind_center + 10]],
                "x",
            )
            plt.xlim(
                (
                    np.min(f1[ind_p[ind_center - 10 : ind_center + 10]]),
                    np.max(f1[ind_p[ind_center - 10 : ind_center + 10]]),
                )
            )

        plt.xlabel("Detuning (MHz)")
        plt.tight_layout()
        plt.show()

    return f


def q_fit_lz(
    f,
    y_res,
    y_mzi,
    fsr_mzi=5.87,
    lw_guess=None,
    wl_res=1550e-9,
    wg_ng=1.50,
    fit_model=0,
    add_drop=False,
    fit_plot=True,
    print_result=True,
):
    """Lorentzian lineshape fitting for a resonance (add-thru / add-drop).

    Args:
        f: freq detuning [MHz]
        y_res: normalized resonance data
        y_mzi: MZI trace
        fsr_mzi: FSR of the MZI fringes in MHz
        lw_guess: starting points for parameter fitting in MHz
        wl_res: resonance wavelength (meters)
        wg_ng: group index
        fit_model (int): 0 non-split, 1 split resonance
        add_drop (bool): treat as add-drop with identical buses
        fit_plot (bool): produce plot
        print_result (bool): print Q fitting results

    Returns:
        q_data (dict): fitted Q, loss, ER, etc.
    """
    if lw_guess is None:
        lw_guess = [3.0, 3.0]

    ind = np.argmin(y_res)
    f = f - f[ind]

    w0 = 0.0
    level = 1.0

    if fit_model == 0:
        fit_str = "add-thru resonator model, no splitting"
        fit_results, _ = curve_fit(TAddThru, f, y_res, p0=[lw_guess[0], lw_guess[1], w0, level])
        lw_fit = fit_results[0:2]
        w0 = fit_results[2]
        level = fit_results[3]
        f = f - w0
        y_res = y_res / level
        y_fit = TAddThru(f, lw_fit[0], lw_fit[1])

    elif fit_model == 1:
        fit_str = "add-thru resonator model, splitting"
        fit_results, _ = curve_fit(
            TAddThruSplit, f, y_res, p0=[lw_guess[0], lw_guess[1], lw_guess[2], w0, level]
        )
        lw_fit = fit_results[0:2]
        w0 = fit_results[3]
        level = fit_results[4]
        f = f - w0
        y_res = y_res / level
        y_fit = TAddThruSplit(f, lw_fit[0], lw_fit[1], fit_results[2])

    if add_drop:
        # Useful when two bus couplings are identical: fitted r_ex is the bus
        # coupling and remaining "intrinsic" loss is r_in + r_ex.
        fit_str = "Add-drop is viewed as add-thru resonator model, no splitting"
        fit_results, _ = curve_fit(TAddThru, f, y_res, p0=[lw_guess[0], lw_guess[1], w0, level])
        lw_fit = fit_results[0:2]
        w0 = fit_results[2]
        level = fit_results[3]
        f = f - w0
        y_res = y_res / level
        y_fit = TAddThru(f, lw_fit[0], lw_fit[1])

    if fit_plot:
        plt.figure()
        plt.plot(f, y_res, label="Resonance")
        plt.plot(f, y_mzi, label=f"{fsr_mzi:.2f} MHz MZI")
        plt.plot(f, y_fit, "--", label="Fit")
        plt.xlim((min(f), max(f)))
        plt.ylim((0, 1.1))
        plt.xlabel("Detuning (MHz)")
        plt.ylabel("Transmission")
        plt.legend()
        plt.title("Loss rates = " + f"{lw_fit[0]:.2f}, {lw_fit[1]:.2f}" + " MHz")
        plt.tight_layout()
        plt.show()

    fwhm = sum(lw_fit)  # MHz; FWHM should not include splitting
    Q_L = c / wl_res / (fwhm * 1e6)
    Q_in = c / wl_res / (lw_fit * 1e6)
    loss_dB = (lw_fit * 2 * np.pi * 1e6) * wg_ng / c * 4.34  # dB/m
    ER = 10 * np.log10(1 / np.min(y_fit))
    g = 0 if fit_model == 0 else fit_results[3]

    if print_result:
        print("-----------Q Analysis------------")
        print("Loss rates = %.2f ," % lw_fit[0], "%.2f MHz" % lw_fit[1])
        print("Total loss rate = %.2f MHz" % fwhm)
        print("Intrinsic Q = %.2f " % (Q_in[0] / 1e6) + "%.2f M" % (Q_in[1] / 1e6))
        print("Loaded Q = %.2f M" % (Q_L / 1e6))
        print("Splitting 2g = %.2f MHz" % (2 * g))
        print("Propagation loss = %.2f," % loss_dB[0] + "%.2f dB/m" % loss_dB[1])
        print("ER = %.2f dB" % ER)

    return {
        "fit info": fit_str,
        "MZIFSR": fsr_mzi,
        "splitting g": g,
        "f detuning": f,
        "y res": y_res,
        "y fit": y_fit,
        "y mzi": y_mzi,
        "loss rates": lw_fit,
        "fwhm": fwhm,
        "intrinsic Q": Q_in,
        "loaded Q": Q_L,
        "ER": ER,
        "loss": loss_dB,
    }
