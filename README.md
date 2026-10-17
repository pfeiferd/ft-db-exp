**Genestrip-FT-DB-Exp** - Experiments concerning
[Genestrip-FT](https://github.com/pfeiferd/genestrip)
===============================================

This project produces the quality and performance results of the Genestrip-FT paper, i.e. the paper
on refining a *k*-mer database's taxonomy so that reads can be classified more specifically than the
original taxonomy allows. Every measurement is taken twice, once on the unrefined database and once
on the refined one, so that any difference is attributable to the refinement alone.

KrakenUniq and Kraken 2 appear in one limited role: given the same genomes as the complete viral
database `cv` and the same read sets, how far does each narrow the species down on the reads it
itself leaves at a genus. That says where another tool's baseline lies, so the refinement's gain
rests on a starting position other tools share.

## License

[Genestrip-FT-DB-Exp is free for any kind of use.](./LICENSE) However, the associated software,
[Genestrip](https://github.com/pfeiferd/genestrip), has a [more restrictive
License](https://github.com/pfeiferd/genestrip#license).

## Building

Needs [Maven 3](https://maven.apache.org/) and a JDK 11 or higher. `mvn install` compiles the
experiment code in `src/main/java`: `org.metagene.ftdbexp` holds the entry points the `exec`
executions call, `…eval` the reports that score a classification, `…nanosim` the genome list NanoSim
simulates from, `…store` the *k*-mer store benchmark. Everything else is driven by the scripts in
`bin` and by Genestrip's own goals, invoked through the `exec` executions of `pom.xml`.

The property `genestrip.version` in `pom.xml` names the version of `org.genestrip:base`, `core` and
`ft` that is resolved, currently **4.0**. A version not on Maven Central has to be put in the local
repository by a `mvn install` in a genestrip checkout. A stale one there is worth ruling out first
when a goal fails on a configuration key it should know: `data/config.properties` sets
`strictConfigCheck=true`, so an unknown key ends the run.

## Requirements

A large disk -- the databases, the RefSeq downloads and the simulated fastq files add up to well
over a terabyte -- and **enough RAM to build the databases**: the peak of one generation run ranges
from 9 GB for `viral` to 46 GB for `nocardia` (`db_ram_mb` of `results/db_gen_perf.csv`), the cgroup
around it up to 57 GB. `pom.xml` caps the heap at `gs.xmx=56G`; the experiments ran on a 64 GB
machine. Matching holds only the *k*-mer store, so the classification experiments are far less
demanding.

The experiments were run against RefSeq Release 233. A later release brings slightly different
results.

## Everything at once

The sections below are the steps in the order they have to run, and that order is also encoded in a
script:

```sh
sh ./bin/run_all_exps.sh              # everything, from the build to the paper's results folder
sh ./bin/run_all_exps.sh --list       # the steps, without running any
sh ./bin/run_all_exps.sh --dry-run    # print the commands instead of running them
sh ./bin/run_all_exps.sh --from 9     # resume at a step after a failure
sh ./bin/run_all_exps.sh --only 4     # a single step
```

Re-running is safe: every step skips work it has already done. Wall time is days. It installs
nothing -- the installation is a prerequisite and needs `sudo`.

| Variable | Effect |
| --- | --- |
| `FRESH_DBS=1` | Run `clean_all.sh` first, so that everything is generated from scratch and the generation is actually measured. Without it an existing database and its `db_gen_*.log` are kept, since a goal that did nothing must not overwrite a real measurement. |
| `FETCH_SALIVA=1` | Fetch the five human saliva runs (about 330 GB) first. Off by default. |
| `KEEP_TICK_SIMS=1` | Never delete a simulated tick fastq. By default those whose abundance table was not preserved are deleted, since only re-simulating produces it -- see "Simulated reads". |
| `SKIP_BUILD=1` | Do not run `mvn install` first. |
| `PAPER_RESULTS` | An existing folder to copy the results to at the end. Unset by default and with no default path: this project does not know where a consumer keeps its inputs. |

Two things it does that are easy to miss by hand. It builds the Java first, because the reports'
labels and columns live in `src/main/java`. And it keeps the order the reports depend on: the
saliva-matched read set before the real-read step, because that is what produces the calibration the
real runs consume, and the Kraken classification before the evaluation that picks its output up.

`FRESH_DBS=1 sh ./bin/run_all_exps.sh` is therefore the from-scratch run: step 3 deletes everything
the experiments ever generated, `results/` included, so no CSV of an earlier batch can survive, and
the last step collects what the paper reads into `paper_results.zip`. The downloads survive --
`clean_all.sh` keeps `data/common` and every fastq a map declares by URL -- so the saliva, tick and
*Nocardia* reads are not fetched again.

## Installation

```sh
sh ./bin/install_tools.sh
MINIMAP2_INDEX_SIZE=24G sh ./bin/install_tools.sh    # faster NanoSim training, more memory
```

Installs cgmemtime (wall time and peak RAM), InSilicoSeq, NanoSim with the binaries it shells out to
(minimap2, LAST, samtools, genometools, bedtools, sam2pairwise), Kraken 2 and KrakenUniq. Everything
lands below `./tools` except a few distribution packages. Idempotent. The real sequencing runs need
nothing from here -- they arrive as gzipped fastq over HTTPS, so sra-toolkit is not installed.

NanoSim departs from the first paper's conda recipe, whose `requirements.txt` pins versions
targeting Python 3.7 that do not build against a current interpreter. Should it break on a numpy or
scikit-learn API change, run it in that conda environment on another machine rather than patching
around it.

## Databases

```sh
sh ./bin/run_exps.sh
```

Downloads the RefSeq genomes, builds the unrefined and the refined database of every project under
`data/projects`, measures the generation performance, computes the intrinsic quality measures --
node precision, subtree precision, restricted subtree precision -- and the tree figures, all into
`./results`.

`strepto`, `tick-borne`, `parasites` and `protozoa` list the human genome in their `additional.txt`,
so that the LCA update accounts for it. Taxid 9606 must **not** also appear in `taxids.txt`, or the
database would store the human genome instead of merely accounting for it. Without the entry every
*k*-mer such a database shares with the human genome stays stored under its own organism and turns
human reads into hits.

`viral` is the exception: its entry is commented out and it sets `maxDust=-1`. Both make `cv`
comparable with the KrakenUniq and Kraken 2 databases over the same genomes, since neither tool can
be given a host genome and neither masks the way the dust filter does. The measured cost is recorded
in the commented block of `data/projects/viral/additional.txt`. For an application rather than a
comparison both belong on.

## Simulated reads

```sh
sh ./bin/make_fastqs.sh                           # both parts
sh ./bin/make_fastqs.sh viral                     # MiSeq (301 bp), HiSeq (126 bp), error-free
ERROR_FREE=1 sh ./bin/make_fastqs.sh protozoa     # the error-free set alone
N_READS=10k sh ./bin/make_fastqs.sh viral         # smoke test; the default is 1M reads per set
```

"MiSeq" and "HiSeq" differ in read length *and* error profile at once, so the error-free set is what
separates the two effects: same length, no errors.

Two further regimes are **not** part of `all` but *are* reported, and the real-read step cannot be
evaluated without them, since a calibration belongs to the pair of database and read set:

```sh
ERROR_SALIVA=1 sh ./bin/make_fastqs.sh viral      # 101 bp at 2.07 %, matched to the saliva runs
ERROR_SALIVA=1 sh ./bin/make_fastqs.sh strepto    # the same, for the second database on those runs
ERROR_MNGS=1 sh ./bin/make_fastqs.sh nocardia     # 70 bp at 2.05 %, matched to the BGISEQ runs
```

No stock InSilicoSeq model reproduces either instrument, so the `basic` model is used at those
parameters. Each regime writes its own fastq map, and each has to be its own invocation: `make_iss`
picks the first non-empty regime variable, so `ERROR_SALIVA=1 … all` would write the saliva-matched
set where the plain Illumina one belongs. `ERROR_NANOPORE=1` (10 %, 125 bp) and
`ERROR_NANOPORE_LONG=1` (10 %, 3,926 bp) also exist and are not reported.

For `tick-borne`, NanoSim trains an error model on the real Nanopore reads of a tick sample and
applies it to the genomes of the twelve tick-borne genera, one simulated fastq per real one for all
eight ticks. Two things are needed first -- the real reads, and the genome list, which labels every
genome `<taxid>x<index>` so that a simulated read's origin stays recoverable from its name:

```sh
mvn exec:exec@fastqdl -Dname=tick-borne -Dfqmap=ticks_real.txt   # tick1..tick8.fastq.gz, ~3.2 GB
mvn exec:exec@nanosimlist -Dname=tick-borne                      # needs the database to exist
```

A fastq map has no checksum column, so `make_fastqs.sh tick-borne` runs `gzip -t` over the downloads
before use and fetches what is missing; `SKIP_FETCH=1` suppresses that. If NCBI's SRA endpoint is
unreachable, `ticks_real.txt` lists two alternatives.

Each tick simulation leaves three files in `results` that cannot be recovered afterwards, since
NanoSim writes every sample into one working directory that is cleared between samples:
`tick-borne_<sample>_quantification.tsv` (the abundance table, derived by aligning the *real* reads
against the genome list), `tick-borne_<sample>_error_rate.tsv` (the per-base error from the
alignment, not from the quality strings -- these runs carry placeholders) and
`tick-borne_nanosim_genomes.tsv`.

**A sample whose `_sim.fastq` is already there is not simulated again**, so no re-running produces
its abundance table. Delete `data/fastq/<sample>_sim.fastq` to force it -- which is what step 5 of
`run_all_exps.sh` does, for exactly those samples whose table is missing.

## Real reads, without ground truth

Three collections: the five human saliva runs, matched against `cv` and against `strepto`; the eight
tick samples, matched against `tb`; and the fourteen clinical metagenomes of ENA study PRJEB34974,
matched against `nocardia`. No ground truth exists for any of them, so what is measured is how far
each database variant narrows the species down.

```sh
sh ./bin/fetch_saliva.sh                  # all five runs, from the ENA over HTTPS, resumable and md5'd
sh ./bin/fetch_saliva.sh SRR5571991       # a single run
mvn exec:exec@fastqdl -Dname=nocardia -Dfqmap=nocardia_mngs.txt    # the Nocardia runs, by ENA URL
sh ./bin/run_classification_exps.sh real  # the evaluation
```

The tick reads are already in place from the simulation above. One mate per saliva run is fetched
and analysed, as in the first paper: 605 to 981 million read pairs each, about 330 GB of gzip for
the five. The first paper used the three `SRR` runs and left the two `ERR` ones of the same study
out; all five are fetched here, so the saliva tables carry five rows. `prefetch`/`fasterq-dump` are
not used -- on these runs they need some 615 GB transiently to arrive at 139 GB and neither step is
resumable.

**One saliva run is measured but not reported.** `SRR5571990` is classified like the rest and stays
in every CSV, but `bin/paper_stats.sh` and `bin/perf_scenarios.sh` leave it out of their averages,
each through an `EXCLUDED_KEYS` set near the top of the Python part, and the paper filters it out of
its tables. Against `strepto` that run leaves 51 reads at a genus where the other four leave one to
thirteen million, which no property of the sample explains; the paper therefore reports four saliva
runs, as the first study on Genestrip reported three. Delete the two sets and the LaTeX filters to
get it back.

The run writes `results/<db>_<key>_specificity.csv`, one row per fastq key: `fastq key`, `sample`,
`reads`, `classified unrefined`, `classified refined`, `genus only`, `genus only share`,
`ungated precision unrefined`, `ungated precision refined`, `rho u`, `rho f`, `est prec g u`, `est
prec g f`.

**It needs a simulated run first.** No ground truth exists to gate with, so the two factors come
from the `rho` columns of a simulated run's `_summary.csv`, named by the last argument of
`run_real`: `iss_saliva` for the saliva runs of `cv` and of `strepto`, `nanosim` for the ticks,
`iss_mngs` for the Nocardia metagenomes. Run `real` before that summary exists and the four
right-hand columns come out empty -- or, worse, are filled from a stale copy still lying in
`results`. A calibration of exactly one row is applied to every sample regardless of its key, which
is how one saliva-like read set stands in for all five runs; more than one row means the rows are
per sample and only a key match will do, which is how each tick borrows the calibration of its own
sample.

The difference between the two ungated precisions is **not** a bound on the gain. It becomes an
estimate of it only through those factors: each `est prec g` column carries one level from what can
be observed to what it stands for.

### The Kraken tools on the real saliva runs

Step 9 of `run_all_exps.sh`, and it has to come before every step that evaluates a read set, i.e.
before steps 10 and 11: the reports pick the per-read output up from `results/kraken` only if it is
already there. It sat after step 10 until 2026-10-17, which silently left the two Kraken rows of the
saliva-like regime at whatever an earlier run had written. The two commands below are
for redoing this part on its own.

```sh
sh ./bin/run_classification_exps.sh kraken     # builds, then classifies simulated AND real
sh ./bin/run_classification_exps.sh real       # evaluates, picking the TSVs up by itself
```

`bin/kraken_build.sh` builds over the very fasta files `cv` is built from. Masking is off
(`--no-masking`), as in the first study on Genestrip, and KrakenUniq is given `--kmer-len 31` to
match Genestrip's *k* while Kraken 2 keeps its defaults -- *k* = 35 over minimizers of 31 with 7
positions spaced out, i.e. a seed of weight 24. `bin/kraken_classify.sh` classifies, with
`FQMAP=<map>` switching it from the four simulated sets to the runs a map names.

**The per-read output is filtered as it is produced**: the tool writes into a FIFO and an `awk`
keeps three fields of the lines whose taxon is not 0. Unfiltered it would be 552 GB per tool.
KrakenUniq would do it itself with `--only-classified-output`, but Kraken 2 parses that option and
then ignores it, so one filter serves both. It costs the measure nothing -- a read left unclassified
is not one left at a genus -- and only `reads` needs repair, which `SpecificityReport` takes from
the Genestrip run.

`SpecificityReport` scores whatever it finds under `results/kraken` and writes
`<db>_<tool>_<key>_specificity.csv` beside its own, each tool on its own subset, the same rule
`RefinementAccuracyReport` follows on the simulated sets. `-Dgs.externalonly=true` on the
`specificity` run redoes those rows alone once the TSVs are in place. The calibration is the tool's
own; its `rho f` is empty, so the refined cells stay empty too.

## Classification quality experiments

```sh
sh ./bin/run_classification_exps.sh
```

Every simulated fastq file is matched against both variants and scored against the reads' known
ground truth, giving one row per fastq file and variant in `results/<db>_<key>_accuracy.csv`: the
boolean positive counts at the genus and species rank, the resolution-based species count -- the
measure that can tell a refined database from an unrefined one at all -- and `unresolved`, the reads
whose ground truth could not be recovered from their name.

Its last six columns serve the ground-truth-free estimate. They are all over one subset
*R<sub>g</sub>*, named by `genus only`: the reads the unrefined database left at or below a genus and
above any species, read off the assignment alone. `genus only score`, `genus only precision species
cand`, `genus only species share` and `genus only zero scoring` are the **gated** measures over it,
which consult the read's true species and so cannot be formed on a real sample; `genus only ungated
precision` is the **ungated** one, the single measure that never consults the ground truth. Both
precisions divide by `genus only`, so their ratio is a statement about the gate alone.

Until 2026-10-13 the subset carried a correctness condition -- the assignment had to be an ancestor
of the read's true species -- and the ungated measure had an observable subset of its own. That
denominator took in the reads a refinement moves off the right path and left out the ones it rescues
from a foreign branch, so it booked the refinement's losses and not its wins.

The calibration the paper transfers to a real sample is a **pair** of factors on the levels, not one
factor on the gain: `rho u` carries the unrefined ungated precision to the gated one it stands for
and `rho f` does the same after the refinement. Both are computed in `_summary.csv`, where the two
variants sit in one row. Ratios of levels rather than of gains, because they then stay within [0, 1]
and because a real sample can be given an estimated precision before and after rather than only an
estimated gain.

### The files the paper's tables read

`_accuracy.csv` is the raw record, one row per fastq file *and variant* -- the wrong shape for a
table, which puts the two variants side by side. Three files are written beside it, one row per read
set:

| File | Holds |
| --- | --- |
| `<db>_<key>_summary.csv` | The restricted precision of both variants over both subsets, and the pair `rho u`, `rho f`. The gains get no column: each is the difference of two columns already there. |
| `<db>_<key>_quality.csv` | The boolean counts at the genus and species rank -- for the *unrefined* variant only, since none of them moves by more than 0.0007 between the two -- and the resolution-based precision, recall and F1 for both, in columns named `… species cand`. |
| `<db>_<key>_simdata.csv` | What the read sets are: read length, per-base error, reads generated, how many fall within the database's scope, how many were unresolved. |

The first three columns of `_simdata.csv` cannot be recovered at evaluation time. `make_fastqs.sh`
measures them as it generates each set, leaves them in `results/<db>_simparams.csv`, and they are
joined in by fastq key. A read set with no row there still gets a row, with those columns empty.

Every file carries the raw `fastq key` beside the printed label, and reports join on the raw key, so
that renaming a label cannot silently break a join.

The performance part also preserves the match goals' own output as
`results/<goal>_<logkey>_<key>.csv`, which says per taxon how many reads were assigned to it. The
log key is load-bearing: the simulated and the real tick runs use the same keys `tick1 .. tick8`
from two different maps. These ride on the performance measurement, so **no cgmemtime means no
per-taxon results**: `run_perf` returns early without it.

### The k = 24 control of cv

A second complete strand for `cv` alone, on the same genomes, taxonomy and read sets, at *k* = 24
instead of 31 -- where Kraken 2's spaced seed of weight 24 sits. It is also the only knob that moves
the share of *k*-mers stored above the data taxa, the mass a refinement can push down, without
touching the genomes, the read sets, the measures or the code. Step 12 of `run_all_exps.sh`.

```sh
sh ./bin/k24_exps.sh all          # the twin, its database, and the four read sets against it
sh ./bin/k24_exps.sh projects     # only the twin project, seconds
sh ./bin/k24_exps.sh build        # only the database and the reports over it
sh ./bin/k24_exps.sh accuracy     # only the classification
K=27 sh ./bin/k24_exps.sh all     # viral-k27 instead
```

`bin/k24_projects.sh` writes the twin: symlinks to the original's `taxids.txt`, `additional.txt` and
`categories.txt`, a `config.properties` from the original's with `kMerSize=24` appended, and the
fastq maps linked under the twin's own name -- which is what makes the two rows comparable, since
the twin is then scored on the very read files `viral` was scored on. `CSV_SUFFIX=-k24` keeps the
timings and disk sizes in CSVs of their own. About 36 minutes and 5.4 GB; it simulates no reads.

### The sampled control of cv

The second twin takes the store down along with the seed weight. Kraken 2 keeps a 4-byte cell per
minimizer and enters only minimizers, a density of 1/3 at its defaults; Genestrip enters every
*k*-mer in a 64-bit word, 50 bits of it in use at *k* = 24, plus a Bloom filter of 10 bits per
entry. Fewer bits per entry and fewer entries, both by about a factor of two, is why a Genestrip
database of Kraken 2's size needs a sampling rate of about four.

```sh
K=24 sh ./bin/sampling_exps.sh all        # the twin, its database, and the four read sets
K=24 sh ./bin/sampling_exps.sh projects   # only the twin project, seconds
```

`cv_k24_s4` comes to 689 MB on disk against the 674 MB of Kraken 2's own hash table and taxonomy
over these genomes, and it is the one variant whose seed weight and whose entry count are both
comparable with Kraken 2's. The sampling selects by the *k*-mer and not by its position, which keeps
the tax ids right: a *k*-mer is entered in every genome it occurs in or in none. A read therefore
keeps about one *k*-mer in four at random positions. `bin/sampling_projects.sh` writes the twin from
the original project, never from the `-k24` twin. Without `K` it builds `viral-s4` at *k* = 31 and
`S=5` sets another rate, but no table reads either.

### The two *k*-mer stores

The appendix compares `RadixKMerStore`, the default, against the `KMerSortedArray` of earlier
versions.

```sh
sh ./bin/store_bench.sh strepto                       # the smallest saliva run, filter on
FILTER=both sh ./bin/store_bench.sh strepto           # with and without the filter
ALL=1 FILTER=both sh ./bin/store_bench.sh tick-borne  # all eight tick runs
```

It needs no second database: it loads the project's refined one, copies it into a sorted array in
the same JVM and classifies the same reads through both, one consumer thread throughout so that the
lookups and not the fastq reader are the limit. It writes `results/storebench_<db>.csv`, one row per
store. `FILTER` takes `on`, `off` or `both`, since the Bloom filter answers most *k*-mers of a real
read set and a run with it on measures something else than one without it.

`perf_scenarios.sh all`, which `run_all_exps.sh` runs, ends with exactly the two calls the appendix
table needs; `perf_scenarios.sh store` runs just those. Mind that the sorted array indexes its
values with a short and so holds 65,535 taxa -- enough for these databases, and a fill beyond that
is refused rather than truncated.

## Machine description

Run **on the machine that executed the experiments**, ideally right after a full run:

```sh
sudo sh ./bin/sysinfo.sh
```

It writes `results/sysinfo.txt`, the full report, and `results/sysinfo.tex`, the facts the paper
quotes as macros. `sudo` only affects the report: DMI, memory-slot and SMART details are root-only,
while the LaTeX file comes from `lscpu`, `/proc`, `lsblk` and `java`. Anything undetermined becomes
`\sysUnknown`, which typesets as a bold **??**, as does an absent file -- so a **??** in that
paragraph means the script has not run on the right machine yet. The JVM reported is the one on
`PATH`, so run it from the same shell as the experiments.

`\sysWorkers` is `\sysThreads - 1`, because `pom.xml` sets `gs.threads` to `-1`, which Genestrip
reads as one thread per processor less one; `AccuracyEvaluator` sets the same value in code, so
every experiment runs with the same thread count. `-Dgs.threads=19` pins it, and the macro then has
to be set by hand.

## Results

Everything is written to `./results`, to be consumed as it is: a paper includes the CSV, SVG and
LaTeX fragments by name from a results folder of its own. `run_all_exps.sh` copies them there as its
last step when `PAPER_RESULTS` names the folder.

**Copying by hand means one more command.** `dbstats.tex` holds every figure the text states as a
macro and is derived from the CSVs rather than written by a goal, and an rsync does not regenerate
it:

```sh
sh ./bin/paper_stats.sh /path/to/the/paper/results
```

Otherwise its tables keep the numbers of the previous batch while the CSVs beside them are new,
which is exactly the mismatch that script exists to prevent. `matchperf.csv` needs nothing extra:
`run_all_exps.sh` writes it into `./results` before the copy.

| File | Written by | Section |
| --- | --- | --- |
| `<db>_dbinfo.csv`, `<db>_ftdbinfo.csv`, `<db>_dbquality.csv`, `<db>_ftquality.csv` | `run_exps.sh` | Databases |
| `<db>_kmerrankstatscsv.csv`, `<db>_branchhistorankcsv.csv` | `run_exps.sh` | Databases |
| `db_gen_<db>.log`, `ftdb_gen_<db>.log` | `run_exps.sh` | Databases |
| `db_gen_perf.csv`, `db_disk_sizes.csv` | `db_gen_perf.sh`, `db_disk_sizes.sh`, from `run_exps.sh` | Databases |
| `<db>_simparams.csv`, `tick-borne_<sample>_quantification.tsv`, `_error_rate.tsv`, `tick-borne_nanosim_genomes.tsv` | `make_fastqs.sh` | Simulated reads |
| `<db>_<key>_accuracy.csv`, `_summary.csv`, `_quality.csv`, `_simdata.csv` | `run_classification_exps.sh` | Classification quality |
| `<db>_<key>_specificity.csv` | `run_classification_exps.sh` | Real reads |
| `<db>_<tool>_<key>_specificity.csv`, `_quality.csv`, `_summary.csv` | `run_classification_exps.sh kraken`, then `real` / `accuracy` | Real reads |
| `match_<logkey>.log`, `ftmatch_<logkey>.log`, `<goal>_<logkey>_<key>.csv` | `run_classification_exps.sh` | Classification quality |
| `perf_<scenario>_<key>_<goal>.log`, `matchperf.csv` | `perf_scenarios.sh` | Classification quality |
| `viral-k24_*`, `db_gen_perf-k24.csv`, `db_disk_sizes-k24.csv` | `k24_exps.sh` | Classification quality |
| `viral-k24-s4_*`, `db_gen_perf-k24-s4.csv`, `db_disk_sizes-k24-s4.csv` | `K=24 sampling_exps.sh` | Classification quality |
| `storebench_<db>.csv`, `storebench_<db>.log` | `store_bench.sh` | Classification quality |
| `dbstats.tex` | `paper_stats.sh` | above |
| `sysinfo.txt`, `sysinfo.tex` | `sysinfo.sh` | Machine description |
| `../paper_results.zip`, beside `results` rather than in it | `collect_paper_results.sh` | below |

**Only the files the paper reads**, as a zip. This is step 16, the last of `run_all_exps.sh`, so a
full run leaves `./paper_results.zip` behind on its own; the commands are for building it again.

```sh
sh ./bin/collect_paper_results.sh                 # ./paper_results.zip, paths kept as results/...
sh ./bin/collect_paper_results.sh --paper-only    # without the paper_stats.sh inputs
OUT=handover.zip sh ./bin/collect_paper_results.sh
```

It is also the closing check of a full run: every name on its list that is not in `results/` is
reported, so a step that silently produced nothing shows up here instead of as an empty table in the
manuscript.

`results` accumulates everything the experiments ever produced and barely a fifth of it reaches the
manuscript. The list of names is fixed inside the script, deliberately: this project does not know
where a paper lives and must not take its path as input. It needs an edit whenever a table gains or
drops a source file -- but nothing goes wrong quietly if that is forgotten, since the script names
every entry it cannot find.

## Cleaning up

```sh
sh ./bin/clean_all.sh                                       # every project, the simulated reads too
sh ./bin/clear_projects.sh viral viral-k24 viral-k24-s4     # one project family
DRY_RUN=1 sh ./bin/clear_projects.sh viral
```

`clean_all.sh` removes the generated databases of every project; the downloaded RefSeq and Genbank
data survives, since re-downloading it takes far longer than rebuilding from it. It is the right
tool when results have to come from one code state, and far too much after a change to a single
project: it takes the simulated reads with it, and new reads would make the rebuilt numbers
incomparable with the tables that were not rebuilt.

`clear_projects.sh` clears `csv/`, `log/`, `krakenout/` and -- what `ftclear` leaves alone -- `db/`,
plus the collected `results/<project>_*`, and keeps the extracted genomes, the simulated fastq files
and `results/kraken`. Deleting is what makes a changed `config.properties` take effect at all: a
Genestrip goal whose output file is in place reports itself made and does nothing.
