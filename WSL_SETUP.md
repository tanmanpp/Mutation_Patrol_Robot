# WSL Setup

Run the web app from WSL so the backend can call Linux bioinformatics tools
directly.

## 1. Install Conda in WSL

Install Miniconda or Miniforge inside WSL. The launcher expects one of these
folders to exist:

```text
~/miniconda3
~/miniforge3
~/anaconda3
~/mambaforge
```

The app launcher creates a conda environment named `mutation_patrol` with:

```text
python, pip, nodejs, samtools, bcftools, minimap2
```

## 2. Enter the Project

If the project is on the Windows `D:` drive:

```bash
cd /mnt/d/YiLun/parasite_AMR/Mutation_Patrol_Robot
```

For best performance with large FASTQ/BAM files, copy the project and data into
the WSL filesystem, for example under `~/projects/Mutation_Patrol_Robot`.

## 3. Start the App from Windows

Double-click:

```text
start_mutation_patrol_robot.bat
```

It will enter WSL, activate/create the conda environment, build the UI, and
start the server. You do not need to move to the project folder first; the bat
file uses its own location as the project path.

Open:

```text
http://localhost:8000
```

## 4. Start the App Manually from WSL

From the project folder:

```bash
bash scripts/run_app_wsl_conda.sh
```

Open:

```text
http://localhost:8000
```

## Development Mode

If you are editing the UI and want Vite hot reload, use two WSL terminals
instead:

```bash
bash scripts/run_backend_wsl.sh
bash scripts/run_frontend_wsl.sh
```

Then open `http://localhost:5173`.

## Notes

- The app stores uploaded files and outputs under `app_data/`.
- Site Query reuses the BAM from a completed analysis run:
  `app_data/results/<run>/work/bam/merged.sorted.bam`, and writes outputs under
  `app_data/results/<run>/site_query/<query>/`.
- The app should run in WSL, not Windows Python, because it calls
  `samtools`, `bcftools`, and `minimap2`.
- The UI can be opened from Windows Chrome/Edge using `http://localhost:8000`.
- If Uvicorn prints a server URL, use `http://localhost:8000` in the browser.
