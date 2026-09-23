# NCOMMS PZT-on-SiN data repository

Data and analysis notebook accompanying:

**Sub-Doppler rubidium atom cooling using a programmable agile integrated PZT-on-SiN resonator**

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

If you use this data or code, please cite the associated Nature Communications paper (citation to be updated upon publication).
