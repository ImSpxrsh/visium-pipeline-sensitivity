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
ssh YOUR_KERBEROS@orcd-login.mit.edu
```

You will need your Kerberos password and Duo. Everything below runs on the
cluster, not your laptop.

## 2. Get the code and set up Python (once)

```bash
cd ~/orcd/pool
git clone https://github.com/ImSpxrsh/visium-pipeline-sensitivity.git
cd visium-pipeline-sensitivity
module load miniforge/26.7.2-0                       # `module spider miniforge` lists versions
conda create -y -p ~/orcd/pool/py312 python=3.12     # the default Python is 3.14
~/orcd/pool/py312/bin/python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt -c constraints.txt   # pinned versions; see below
mkdir -p logs/slurm
```

`constraints.txt` pins the package versions the test suite passes with. With
the newest pandas and statsmodels, the variance decomposition silently returns
NaN, so don't skip it. Run `python -m pytest -m "not network" -q` to check:
all tests should pass.

If the repo is private, `git clone` asks for a GitHub token, not your
password. Alternatively, copy the repo from your laptop with `scp` or `rsync`.

## 3. Download the data (once per dataset)

Compute nodes on the regular partition may not be able to download, and a job
that tries will hang until it hits its time limit. The mouse-brain loader also
fetches a 398 MB H&E image the first time a pipeline runs, which loading the
dataset alone does not trigger. The simplest way to fetch everything is to run
one spec per dataset on the 15-minute test partition, which can download:

```bash
source .venv/bin/activate
python scripts/make_factorial_config.py
for d in visium_hne breast_cancer; do
  srun -p mit_quicktest --time=00:14:00 --cpus-per-task=4 --mem=16G \
    python scripts/main.py --config configs/multiverse_full.yaml --dataset $d \
    --results-root results/probe_$d --spec-index 0
done
```

Files land in `data/` inside the repo, and later jobs reuse them.

## 4. Test one batch on the regular partition

Each array task runs 8 consecutive specs, so 864 specs is tasks 0-107.

```bash
sbatch --array=0 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
squeue -u $USER                       # watch it; empty means it finished
cat logs/slurm/mosaic-spec_*_0.out    # check for errors
```

The job scripts use the `mit_normal` partition (12-hour limit). If Slurm
rejects it, run `sinfo` and change the `--partition` line in the `.sbatch`
files.

## 5. Submit the full sweep

```bash
sbatch --array=0-107 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
```

`mit_normal` allows at most 448 queued jobs and 96 CPUs per user, which is why
specs are batched 8 per task: 108 tasks per dataset, about 24 running at once.
On Engaging a spec took a median of about 2 minutes (5 minutes at most), and a
whole 864-spec sweep finished overnight.

Slurm prints a job ID. Queue the results step (S3, ARI, headline instability)
and then the analysis step (variance decomposition, stakes, confound, figures)
to run after it:

```bash
sbatch --dependency=afterany:ARRAY_ID scripts/slurm/aggregate.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
sbatch --dependency=afterok:AGG_ID scripts/slurm/analyze.sbatch visium_hne results/full_hne
```

Both use `--lean`, which keeps only cluster labels and SVG tables in memory.
Without it, loading 864 full results (about 150 MB each for mouse brain) runs
out of memory. Add `--mail-type=END,FAIL --mail-user=YOU@mit.edu` to either
`sbatch` command to get an email when it finishes.

Repeat with `breast_cancer` and `results/full_breast_cancer` for the tumor
section.

Each spec is cached by its content hash, so if some jobs fail or time out you
can resubmit just those tasks (`--array=8,22,...`): finished specs are
skipped. Failures for each spec are recorded in
`results/full_hne/runs/failure_manifest_specN.json`.

**Stuck nodes.** In our run, 16 of 216 tasks sat for 12 hours without
finishing one spec and left no log. Most were on `node3105`, `node3109`, and
`node1618`. If a task produces nothing within an hour, resubmit it with a
shorter limit and those nodes excluded:

```bash
sbatch --time=01:00:00 --exclude=node3105,node3109,node1618 --array=8,22 scripts/slurm/sweep_array.sbatch configs/multiverse_full.yaml visium_hne results/full_hne
```

## 6. Storage

Per-spec results take about 150 MB each for mouse brain and 45 MB for breast
cancer, so one full sweep is about 127 GB and 38 GB. On our account the
`pool` volume held only 177 GB, even though the login quota report listed
1 TB, and once it fills up, writes fail. Check with `df -h ~/orcd/pool`. Keep
the code in `pool`, and put per-spec runs in scratch
(`~/orcd/scratch`, 1 TB) with a link back into the repo:

```bash
mkdir -p ~/orcd/scratch/mosaic_runs/full_hne
ln -s ~/orcd/scratch/mosaic_runs/full_hne results/full_hne/runs   # before submitting
```

Scratch is not backed up, so copy anything you need to keep.

## 7. Bring results back

Copy only the aggregates to your laptop:

```bash
# run this on your laptop
scp -r YOUR_KERBEROS@orcd-login.mit.edu:~/orcd/pool/visium-pipeline-sensitivity/results/full_hne/aggregate results/full_hne/
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
