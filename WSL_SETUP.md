# WSL Setup

Mutation Patrol Robot should be run from **Linux**, **Windows with WSL**, or
**macOS with conda/bioconda**. Native Windows Python is not recommended because
the backend calls Linux bioinformatics tools such as `samtools`, `bcftools`,
and `minimap2`.

This guide is for Windows users running the app through WSL.

## Quick Start

1. Install WSL and a Linux distribution, such as Ubuntu.
2. Install Miniconda or Miniforge inside WSL.
3. Double-click this file from Windows:

```text
start_mutation_patrol_robot.bat
```

The launcher will:

- enter WSL automatically
- check that WSL and conda are available
- create or synchronize the conda environment `mutation_patrol`
- verify Python, Node.js, `samtools`, `bcftools`, and `minimap2`
- install backend/frontend dependencies only when needed
- build the web UI only when its source changed
- start the server

When the server is ready, open:

```text
http://localhost:8000
```

You do **not** need to manually move into the project folder before using the
`.bat` file. It detects its own project path.

## Install WSL

Open PowerShell as Administrator and run:

```powershell
wsl --install
```

Restart Windows if prompted. Then open Ubuntu from the Start menu once to finish
creating your Linux user account.

To check that WSL is available:

```powershell
wsl --list --verbose
```

## Install Conda Inside WSL

Install Miniconda or Miniforge inside WSL. The launcher looks for conda in one
of these locations:

```text
~/miniconda3
~/miniforge3
~/anaconda3
~/mambaforge
```

Miniforge is recommended because it works well with `conda-forge` and
`bioconda`.

After installing conda, close and reopen your WSL terminal, then check:

```bash
conda --version
```

## Conda Environment

The environment definition is stored in:

```text
environment.yml
```

The app launcher creates or synchronizes this environment automatically:

```text
mutation_patrol
```

It installs:

```text
python=3.11
pip
nodejs
samtools
bcftools
minimap2
```

You usually do not need to create this environment yourself. The first launch
may take several minutes because packages and frontend dependencies are
installed. Later launches reuse the verified environment.

To verify the active environment manually:

```bash
conda activate mutation_patrol
python scripts/check_environment.py
```

The command exits with an error and lists every missing component when the
environment is incomplete.

## Manual Start from WSL

If you prefer to start the app manually, open WSL and run:

```bash
cd /mnt/d/YiLun/parasite_AMR/Mutation_Patrol_Robot
bash scripts/run_app_wsl_conda.sh
```

Then open:

```text
http://localhost:8000
```

If your project is in a different folder, replace the `cd` path with your own
project path.

## Development Mode

Use this only if you are editing the backend or frontend and want live reload.
Both development launchers run the same conda preflight as the normal app
launcher before starting.

Open one WSL terminal for the backend:

```bash
cd /mnt/d/YiLun/parasite_AMR/Mutation_Patrol_Robot
bash scripts/run_backend_wsl.sh
```

Open another WSL terminal for the frontend:

```bash
cd /mnt/d/YiLun/parasite_AMR/Mutation_Patrol_Robot
bash scripts/run_frontend_wsl.sh
```

Then open:

```text
http://localhost:5173
```

## Where Files Are Stored

Uploaded files and results are stored under:

```text
app_data/
```

Important output folders:

```text
app_data/databases/
app_data/results/
app_data/uploads/
```

Site Query reuses the BAM from a completed sample analysis:

```text
app_data/results/<run>/work/bam/merged.sorted.bam
```

Site Query outputs are written to:

```text
app_data/results/<run>/site_query/<gene>/site_query.csv
```

## Performance Note

For small tests, running from the Windows drive is fine:

```text
/mnt/d/YiLun/parasite_AMR/Mutation_Patrol_Robot
```

For large FASTQ/BAM files, WSL is faster if the project and data are copied into
the Linux filesystem, for example:

```text
~/projects/Mutation_Patrol_Robot
```

## Troubleshooting

### Browser Cannot Open the App

Use this URL in Windows Chrome or Edge:

```text
http://localhost:8000
```

Do not use `http://0.0.0.0:8000` in the browser.

### Conda Is Not Found

Make sure conda is installed inside WSL, not only on Windows:

```bash
conda --version
```

If this fails, install Miniconda or Miniforge inside WSL and reopen the WSL
terminal.

### npm Package Error

The conda package `nodejs` should provide `npm`. The launcher checks this and
will try to install `nodejs` again if `npm` is missing.

### Port 8000 Is Already in Use

Stop the existing server process, or change the port in:

```text
scripts/run_app_wsl_conda.sh
```

### samtools/bcftools/minimap2 Missing

Activate the environment and check:

```bash
conda activate mutation_patrol
samtools --version
bcftools --version
minimap2 --version
```

If any command is missing, reinstall the environment packages:

```bash
conda env update -n mutation_patrol -f environment.yml --prune
python -m pip install -r backend/requirements.txt
```

You can also delete `app_data/.environment.sha256` and start the app again to
force the launcher to synchronize the environment.
