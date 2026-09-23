# NCOMMS PZT-on-SiN data repository

Data and analysis notebook accompanying:

**[Sub-Doppler rubidium atom cooling using a programmable agile integrated PZT-on-SiN resonator](https://doi.org/10.1038/s41467-026-77526-9)**

## Contents

| Path | Description |
|------|-------------|
| `780_PZT_NCOMMS_DATA.ipynb` | Main analysis notebook (Figures 2, 4, 5) |
| `resonators.py` | Resonator Q-fitting helpers used by the notebook |
| `plot_style.mplstyle` | Matplotlib style sheet |
| `figure2/` | Resonator / PZT tuning measurement CSVs |
| `figure4/` | Time-of-flight fluorescence images (TIFF) |
| `figure5/` | Beat-note and MOT atom-number CSVs |
| `pyproject.toml` / `uv.lock` | Python dependencies (managed with [uv](https://docs.astral.sh/uv/)) |

## Setup

Requires [uv](https://docs.astral.sh/uv/) and Python 3.12+.

```bash
git clone https://github.com/aisichenko/ncomms2026-780pzt.git
cd ncomms2026-780pzt
uv sync
```

Register a Jupyter kernel (once):

```bash
uv run python -m ipykernel install --user \
  --name=ncomms-pzt-data \
  --display-name="Python (ncomms-pzt-data)"
```

Open `780_PZT_NCOMMS_DATA.ipynb`, select the **Python (ncomms-pzt-data)** kernel, and run with this folder as the working directory.

Alternatively:

```bash
uv run jupyter lab
```

## Citation

If you use this data or code, please cite:

Isichenko, A., Carpenter, S., Montifiore, N. *et al.* Sub-Doppler rubidium atom cooling using a programmable agile integrated PZT-on-SiN resonator. *Nat Commun* (2026). https://doi.org/10.1038/s41467-026-77526-9

```bibtex
@article{isichenko2026subdoppler,
  title   = {Sub-Doppler rubidium atom cooling using a programmable agile integrated {PZT}-on-{SiN} resonator},
  author  = {Isichenko, Andrei and Carpenter, Steven and Montifiore, Nick and Wang, Jiawei and Dangi, Mayand and Chauhan, Nitesh and Mukherjee, Pritha and Yang, Xuting and Indukuri, Nitin and Harrington, Mark W. and Zhong, Chuan and Kierzewski, Iain M. and Rudy, Ryan Q. and Choy, Jennifer T. and Blumenthal, Daniel J.},
  journal = {Nature Communications},
  year    = {2026},
  doi     = {10.1038/s41467-026-77526-9},
  url     = {https://www.nature.com/articles/s41467-026-77526-9}
}
```
