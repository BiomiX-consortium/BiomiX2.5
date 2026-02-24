# BiomiX – Linux Installation Guide

## Prerequisites

### 1. Conda / Mamba

BiomiX requires a working **conda** installation with **mamba** available as a
solver. If you already have Miniforge or Mambaforge installed, skip this step.

Install [Miniforge](https://github.com/conda-forge/miniforge) (recommended —
ships with mamba pre-installed):

```bash
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
```

Accept the licence, choose an install prefix (default: `~/miniforge3`), and let
the installer initialise your shell. **Restart your terminal** before
continuing.

> **Non-default conda prefix?**
> If conda is installed under a path other than `/opt/conda` or `~/miniforge3`,
> export `CONDA_BASE` before running the install script:
> ```bash
> export CONDA_BASE=/path/to/your/conda
> ```

---

## Installation

### 2. Download BiomiX

```bash
git clone https://github.com/IxI-97/BiomiX.git
cd BiomiX
```

Or download and unpack the release archive manually from
<https://github.com/IxI-97/BiomiX/releases>.

### 3. Run the install script

```bash
bash _INSTALL/BiomiX_INSTALL_LINUX.sh
```

What the script does:

| Step | Action |
|------|--------|
| 1 | Removes any pre-existing `biomix` conda environment |
| 2 | Runs `mamba env create -f _INSTALL/biomix_env.yml` — one solver pass that installs Python, R, all CRAN binaries, all Bioconductor binaries, and PyPI packages |
| 3 | Runs `Rscript _INSTALL/MODULE_LINUX_custom_r.r` — downloads `Package_linux.tar.gz` from GitHub and installs the 8 packages not available on conda (masstools, massdataset, metid, metpath, litsearchr, cmmr, NEMO, MintTea) |

> **Internet access is required** during installation to pull conda packages and
> the `Package_linux.tar.gz` tarball (~500 MB).

The full install typically takes **5–15 minutes** on a modern server.

---

## Launching BiomiX

Activate the environment and start the application:

```bash
source "$CONDA_BASE/etc/profile.d/conda.sh"   # adjust path if needed
conda activate biomix
bash LaunchApp_BiomiX_Linux.sh
```

Or in a single command:

```bash
CONDA_BASE=/opt/conda bash LaunchApp_BiomiX_Linux.sh
```

BiomiX opens a browser-based interface (Shiny). If running on a **headless
server**, forward port 3838 via SSH and open `http://localhost:3838` locally:

```bash
ssh -L 3838:localhost:3838 user@your-server
```

---

## Updating an existing environment

If you already have `biomix` installed and want to sync it with a newer
`biomix_env.yml`:

```bash
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate biomix
mamba env update -n biomix -f _INSTALL/biomix_env.yml --prune
Rscript _INSTALL/MODULE_LINUX_custom_r.r
```

---

## Verifying the installation

```bash
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate biomix

Rscript -e "library('ggpubr'); library('DESeq2'); library('MOFA2'); cat('R OK\n')"
python -c "import pandas; import sklearn; print('Python OK')"
```

Both lines should print `OK` without errors.

---

## Troubleshooting

### `mamba: command not found`

Mamba is bundled with Miniforge. If you used plain Miniconda, install mamba
first:

```bash
conda install -n base -c conda-forge mamba
```

### `CONDA_BASE` not found / conda not initialised

Re-run:

```bash
source /path/to/conda/etc/profile.d/conda.sh
```

or add it to your `~/.bashrc` / `~/.zshrc`.

### MintTea install warning

MintTea is installed from GitHub. If the server has no outbound internet access,
the script will print a warning and continue. All other packages install
normally. MintTea can be installed later manually:

```bash
conda activate biomix
Rscript -e "devtools::install_github('efratmuller/MintTea', upgrade='never')"
```

### Disk space

The full `biomix` environment requires approximately **8–12 GB** of disk space.
