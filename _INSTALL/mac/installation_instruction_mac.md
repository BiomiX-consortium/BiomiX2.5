# BiomiX – macOS Installation Guide

## Requirements

- macOS 11.0 Big Sur or later (Apple Silicon or Intel)
- ~12 GB free disk space
- Internet connection during installation
- No prerequisites for end users — Miniforge is installed automatically if needed

---

## Installation

### 1. Download BiomiX2.5.pkg

Download `BiomiX2.5.pkg` from the [BiomiX releases page](https://github.com/IxI-97/BiomiX/releases).

### 2. Open the installer

Double-click `BiomiX2.5.pkg` to launch the macOS Installer.

> **Gatekeeper warning?**
> If macOS says the package cannot be opened because it is from an unidentified
> developer, right-click (or Control-click) the file and choose **Open**, then
> confirm in the dialog. This is expected for unsigned packages distributed
> outside the App Store.

### 3. Follow the installer steps

Click through **Introduction → License → Installation Type → Install**.
Enter your macOS password when prompted.

### 4. Wait for installation to complete

After the installer copies files, a **postinstall script** runs in the background:

| Step | Action |
|------|--------|
| 1 | Detects or installs **Miniforge3** (if no conda is found) |
| 2 | Installs **mamba** into the base environment (if not present) |
| 3 | Runs `mamba env create -f biomix_env.yml` — one solver pass for Python, R, all CRAN/Bioconductor binaries, and PyPI packages |
| 4 | Runs three R scripts to install packages not available on conda (masstools, massdataset, metid, metpath, litsearchr, cmmr, SNFtool, NEMO, MintTea) |

> **This takes 15–30 minutes.** The macOS Installer may appear to be stuck on
> "Running package scripts" — this is normal. Do not cancel.

When the installer closes, **`BiomiX.app`** is available in `/Applications`.

---

## Launching BiomiX

Double-click **BiomiX.app** in `/Applications` (or Launchpad).

What happens:
1. A **Terminal window** opens automatically.
2. The `Biomix2.5` conda environment is activated.
3. BiomiX starts and your default browser opens at `http://localhost:3838`.

To quit, close the Terminal window that was opened by the app.

---

## Verifying the Installation

Open **Terminal** and run:

```bash
# Check the conda environment exists
conda env list | grep Biomix2.5

# Check R packages
conda run -n Biomix2.5 Rscript -e "library('DESeq2'); library('MOFA2'); cat('R OK\n')"

# Check Python packages
conda run -n Biomix2.5 python -c "import pandas; import sklearn; print('Python OK')"
```

Both commands should print `OK` without errors.

To review the install log for errors:

```bash
grep -i biomix /var/log/install.log
```

---

## Developer: Building the .pkg

### Prerequisites

- macOS with Xcode Command Line Tools:
  ```bash
  xcode-select --install
  ```
- The full BiomiX repository (cloned locally):
  ```bash
  git clone https://github.com/IxI-97/BiomiX.git
  cd BiomiX
  ```

### Build

Run the build script from the **repository root**:

```bash
bash _INSTALL/mac/build_pkg.sh
```

What the script does:

| Step | Action |
|------|--------|
| 1 | Copies BiomiX source files into the app payload (excluding Windows-only files, `.git`, etc.) |
| 2 | Copies `biomix_env.yml` and R install scripts into `BiomiX.app/Contents/Resources/` |
| 3 | Runs `pkgbuild` to produce `BiomiX2.5.pkg` in the repository root |

Output: **`BiomiX2.5.pkg`** in the repository root.

### Testing locally (no signing required)

```bash
sudo installer -pkg BiomiX2.5.pkg -target /
```

### Signing and notarization

For **lab/internal distribution**, signing is not required.
For **public distribution** (e.g. posted on a website), Apple requires the package
to be signed and notarized. See Apple's documentation on
[Notarizing macOS software before distribution](https://developer.apple.com/documentation/security/notarizing-macos-software-before-distribution).

---

## Troubleshooting

### "BiomiX.app is damaged and can't be opened" / Gatekeeper block

Remove the quarantine attribute:

```bash
xattr -cr /Applications/BiomiX.app
```

Then try opening the app again.

### BiomiX.app shows an alert: "conda not found"

The app cannot locate a conda installation in the usual places. Check that
Miniforge was installed during the postinstall step:

```bash
grep -i miniforge /var/log/install.log
ls ~/miniforge3/etc/profile.d/conda.sh
```

If missing, install Miniforge manually:

```bash
# Apple Silicon
curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh | bash
# Intel
curl -fsSL https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-x86_64.sh | bash
```

Then re-run the postinstall script manually:

```bash
sudo bash /Applications/BiomiX.app/Contents/Resources/../../../_scripts_/postinstall
```

Or recreate the environment directly:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
mamba env create -f /Applications/BiomiX.app/Contents/Resources/biomix_env.yml
```

### Install appears frozen

Normal — conda environment creation can take **15–30 minutes**. Do not cancel
the installer. You can monitor progress by watching the install log in a
separate Terminal:

```bash
tail -f /var/log/install.log
```

### Disk space

The full `Biomix2.5` conda environment requires approximately **10–12 GB**
of disk space. Ensure you have at least this much free before installing.

### mamba not found

The postinstall script falls back to `conda env create` automatically if mamba
is not available. The install will still succeed, though it may be slower.

To install mamba manually afterwards:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda install -n base -c conda-forge mamba
```

### MintTea install warning

MintTea is installed from GitHub during postinstall. If the machine had no
outbound internet access at install time, the warning is printed and the
remaining packages install normally. Install MintTea manually later:

```bash
source ~/miniforge3/etc/profile.d/conda.sh
conda activate Biomix2.5
Rscript -e "devtools::install_github('efratmuller/MintTea', upgrade='never')"
```
