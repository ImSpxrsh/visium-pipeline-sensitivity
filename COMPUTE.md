# Running sweeps on the MIT cluster

The 30-spec mouse-brain sweep takes about 3.5 hours on a laptop (one spec took
95 minutes), and the breast-cancer sweep takes longer. The full factorial in
`configs/multiverse_full.yaml` has 864 specs, which would take days in serial.
A cluster runs many specs at once, so the same work finishes in a few hours.

These steps assume a Slurm cluster, such as MIT ORCD's Engaging cluster. Login
hostnames and partition names change, so check the ORCD docs
(https://orcd-docs.mit.edu) or the welcome email you got with your account.

## 1. Log in

```bash
ssh YOUR_KERBEROS@orcd-login001.mit.edu
```

You will need your Kerberos password and Duo. Everything below runs on the
cluster, not your laptop.

## 2. Get the code and set up Python (once)

```bash
cd ~
git clone https://github.com/ImSpxrsh/MOSAIC-ISEF.git
cd MOSAIC-ISEF
module avail python          # list the Python versions the cluster provides
module load python/3.12      # pick a 3.12 (or 3.11) from that list; exact name varies
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
mkdir -p logs/slurm
```

If `module load` fails, try `module load miniforge` and then create the venv
the same way.

## 3. Download the data on the login node (once per dataset)

Compute nodes often have no internet access, so download first:

```bash
source .venv/bin/activate
python -c "from mosaic.data_loader import load_dataset; load_dataset('visium_hne')"
python -c "from mosaic.data_loader import load_dataset; load_dataset('breast_cancer')"
```

## 4. Generate the full-factorial config

```bash
python scripts/make_factorial_config.py      # writes configs/multiverse_full.yaml (864 specs)
```

## 5. Test one spec before submitting 864

```bash
sbatch --array=0 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
squeue -u $USER                       # watch it; empty means it finished
cat logs/slurm/mosaic-spec_*_0.out    # check for errors
```

If Slurm complains about the partition, run `sinfo` to list partitions, then
uncomment the `--partition` line in both `.sbatch` files and set it.

## 6. Submit the full sweep

```bash
sbatch --array=0-863%50 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
```

`%50` limits it to 50 jobs at a time so it doesn't exhaust your allocation.
Slurm prints a job ID. Then queue the aggregation step to run when the array
finishes:

```bash
sbatch --dependency=afterany:JOB_ID scripts/slurm/aggregate.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
```

Repeat with `breast_cancer` and `results/full_breast` for the tumor section.

Each spec is cached by its content hash, so if some jobs fail or time out you
can resubmit the same array: finished specs are skipped. Failures for each
spec are recorded in `results/full_hne/runs/failure_manifest_specN.json`.

## 7. Bring results back

The per-spec `.h5ad` files are large; keep them on the cluster (use scratch
storage if your home quota is small). Copy only the aggregates to your laptop:

```bash
# run this on your laptop
scp -r YOUR_KERBEROS@orcd-login001.mit.edu:~/MOSAIC-ISEF/results/full_hne/aggregate results/full_hne/
```

## What to run, in order of value to the paper

1. **Full factorial on mouse brain and breast cancer.** The paper's variance
   shares come from a curated set centered on one pipeline, which the paper now
   lists as a limitation. A balanced 864-spec run tests whether resolution
   still dominates.
2. **More tissue sections.** The tumor-margin result rests on one breast-cancer
   section. Repeating it on more annotated sections (e.g. the 12 DLPFC
   sections from spatialLIBD, or other tumor sections with pathology labels)
   is the biggest gain in credibility. The second section from the same
   breast-cancer block already has a loader (`--dataset breast_cancer_2`), so
   the boundary result can be checked there right away (the SEDR pathology
   labels only cover section 1). DLPFC and other datasets need new loaders in
   `mosaic/data_loader.py` first.
3. **Vary the SVG test.** Geary's C is already implemented (`gearys_c`) but
   was not in any sweep. Adding it as an axis doubles the factorial to 1728
   specs and lets the paper report how much the SVG method matters.
