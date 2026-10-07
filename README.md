**Genestrip-FT-DB-Exp** - Experiments concerning [Genestrip-FT](https://github.com/pfeiferd/genestrip)
===============================================

This project runs the quality and performance experiments behind the Genestrip-FT paper, i.e. the
paper on refining a *k*-mer database's taxonomy so that reads can be classified more specifically
than the original taxonomy allows.

It is the successor of `genestrip-db-exp`, which produced the results of the first Genestrip paper.
Where that project compared Genestrip against KrakenUniq, Kraken 2 and Ganon, this one compares
Genestrip against *itself*: every measurement is taken twice, once on the unrefined database and
once on the refined one, so that any difference is attributable to the refinement alone. No foreign
classifier is therefore needed here.

## License

[Genestrip-FT-DB-Exp is free for any kind of use.](./LICENSE)
However, the associated software, [Genestrip](https://github.com/pfeiferd/genestrip), has a
[more restrictive License](https://github.com/pfeiferd/genestrip#license).

## Building

The project requires [Maven 3](https://maven.apache.org/) and a JDK 11 or higher. Given a matching
installation, `mvn install` compiles the Java part in `src/main/java`, which holds the experiment
code that is not part of Genestrip itself:

| Package | Purpose |
| --- | --- |
| `org.metagene.ftdbexp.eval` | scores read classifications against the ground truth of simulated reads |
| `org.metagene.ftdbexp.nanosim` | writes the genome list NanoSim simulates from |

Everything else is driven by the shell scripts in `bin` and by Genestrip's own goals, invoked
through the `exec` executions declared in `pom.xml`.

## Requirements

You need a large disk -- the databases, the RefSeq downloads and the simulated fastq files add up to
well over a terabyte -- and, more importantly, **enough RAM to build the databases**. Building the
`viral` database allocates roughly 18 GB in a single block and the larger databases need more; the
original experiments used a machine with 64 GB. A machine that cannot build a database can still
*match* against one, since matching only holds the *k*-mer store in memory, so the classification
experiments are far less demanding than the database generation.

The experiments were run against RefSeq Release 233. A later release brings slightly different
results.

## 0. Everything at once

The sections below are the steps in the order they have to run. That order is also encoded in a
script, which drives them rather than replacing them:

```sh
sh ./bin/run_all_exps.sh              # everything, from the build to the paper's results folder
sh ./bin/run_all_exps.sh --list       # the steps, without running any
sh ./bin/run_all_exps.sh --dry-run    # print the commands instead of running them
sh ./bin/run_all_exps.sh --from 9     # resume at a step after a failure
sh ./bin/run_all_exps.sh --only 4     # a single step
```

Re-running is safe: every underlying step skips work it has already done. Wall time is days, so
`--from` is what a run interrupted on day two is resumed with.

It does **not** install anything -- section 1 is a prerequisite, not a step, since it needs `sudo`
for a few distribution packages. Run `install_tools.sh` once beforehand.

| Variable | Effect |
| --- | --- |
| `FRESH_DBS=1` | Delete the databases first, so that their generation is actually measured. Without it an existing database is kept and its `db_gen_*.log` / `ftdb_gen_*.log` are left alone -- `run_exps.sh` will not overwrite a real measurement with the timing of a goal that did nothing. Set it whenever the performance table is to be rebuilt from one consistent batch. |
| `FETCH_SALIVA=1` | Fetch the human saliva runs (about 400 GB) before the real-read step. Off by default. |
| `KEEP_TICK_SIMS=1` | Never delete a simulated tick fastq. By default the script deletes exactly those whose NanoSim abundance table was not preserved, since only re-simulating produces it -- see section 3. |
| `SKIP_BUILD=1` | Do not run `mvn install` first. |
| `PAPER_RESULTS` | An existing folder to copy the results to at the end, in addition to `./results`. Unset by default, and there is no default path: this project does not know where a consumer keeps its inputs. |

Two things the script does that are easy to miss when running the steps by hand. It builds the Java
first, because the evaluation reports -- the labels their CSV files carry and the columns they are
made of -- live in `src/main/java`, so a run against a stale build reproduces whatever those classes
did last time. And it runs the saliva-matched read set *before* the real-read step, because that is
what produces the calibration the real saliva runs consume; in the other order the estimate columns
come from a stale copy or from nothing. Its preflight also reports a missing `cgmemtime` or fastq map
up front rather than however many hours later the step that needs it would.

## 1. Installation

Everything external is installed below `./tools` and nothing but a few distribution packages is
touched outside the project:

```sh
sh ./bin/install_tools.sh
```

This installs cgmemtime (wall time and peak RAM measurement), InSilicoSeq (Illumina read simulation)
and NanoSim (Nanopore read simulation), including the binaries NanoSim shells out to -- minimap2,
LAST, samtools, genometools, bedtools and sam2pairwise. The script is idempotent, so it can be
re-run after a partial failure.

Two remarks on NanoSim, which is the awkward one:

* The installation deliberately departs from the conda recipe used for the first paper. NanoSim's
  `requirements.txt` pins versions that target Python 3.7 and do not build against a current
  interpreter, and the conda hosts are not reachable from every network. The external binaries
  therefore come from the distribution and the Python dependencies from PyPI without the pins. If
  NanoSim ever breaks on an API change in a newer numpy or scikit-learn, run it in the original
  conda environment on another machine rather than patching around it.
* With a reference larger than minimap2's default index batch size of 4G the index gets split, which
  slows NanoSim's training down considerably. The first paper's experiments therefore raised it. The
  patch is available but **off by default**, because minimap2 then holds the whole index in memory
  and a value beyond the machine's RAM turns a slow run into a failing one:

  ```sh
  MINIMAP2_INDEX_SIZE=24G sh ./bin/install_tools.sh
  ```

## 2. Databases

```sh
sh ./bin/run_exps.sh
```

This downloads the RefSeq genomes, builds the unrefined and the refined database for every project
under `data/projects`, measures the generation performance, and computes the intrinsic quality
measures -- node precision and subtree precision -- along with the tree figures. All resulting CSV,
SVG and LaTeX files are collected in `./results`.

Several projects list the human genome in their `additional.txt`, `viral` among them since `cv` is
applied to human saliva. The taxid there tells the LCA update whose genome the file is; it must
**not** also appear in `taxids.txt`, or the database would store the human genome rather than merely
account for it. Without that entry every *k*-mer a virus shares with the human genome -- endogenous
retroviral sequence, integrated herpesvirus, host contamination in the viral assemblies -- stays
stored under its virus and turns human reads into viral hits, which on a saliva sample is most of
the data. The file itself is declared once in `data/common/fasta/downloads.txt` and shared.

## 3. Simulated reads

```sh
sh ./bin/make_fastqs.sh              # both parts
sh ./bin/make_fastqs.sh viral        # InSilicoSeq only
sh ./bin/make_fastqs.sh tick-borne   # NanoSim only
```

For `viral`, InSilicoSeq applies its Illumina "MiSeq" and "HiSeq" error models to all RefSeq genomes
of the category "Viral", one million reads each as in the first paper. `N_READS=10k` produces a much
smaller set for a quick smoke test.

`all` generates the two Illumina sets and the error-free one of every InSilicoSeq project, which is
what the paper reports. "MiSeq" and "HiSeq" differ in read length *and* error profile at once, so the
error-free set is what separates the two effects: it keeps the read length and drops the errors.

```sh
sh ./bin/make_fastqs.sh protozoa                          # MiSeq (301 bp) and HiSeq (126 bp)
ERROR_FREE=1 sh ./bin/make_fastqs.sh protozoa             # "perfect": same length, no errors
```

A fourth regime is **not** part of `all` but *is* reported, and the real saliva runs of step 4 cannot
be evaluated without it:

```sh
ERROR_SALIVA=1 sh ./bin/make_fastqs.sh viral              # 101 bp at 2.07 % per-base error
ERROR_SALIVA=1 sh ./bin/make_fastqs.sh strepto            # the same, for the second database
```

Two databases are classified against the saliva runs -- `cv` and `strepto`, saliva being a
streptococcal habitat -- and a calibration belongs to the pair of database and read set, so each
needs a saliva-matched set of its own. The accuracy runs over both are what
`sh ./bin/run_classification_exps.sh saliva` does, in one call and with the regime variable set and
restored inside it.

The real saliva runs are Illumina HiSeq 2000 at 101 base pairs and 2.07 % per-base error, which no
stock InSilicoSeq model reproduces -- its "HiSeq" model is an order of magnitude cleaner than the
instrument it is named after. This set is generated by the `basic` model at those parameters instead,
so that the calibration factors of the paper's Section "Estimating the gain without ground truth" can
be *measured* at the parameters of the data they are applied to rather than extrapolated to them.

It has to be its own invocation. `make_iss` picks its regime from the first regime variable that is
non-empty, so `ERROR_SALIVA=1 sh ./bin/make_fastqs.sh all` would write the saliva-matched set where
the plain Illumina one belongs -- under a different mapping file, so the run would look like it had
succeeded. Both `*_all_regimes` functions therefore clear every regime variable before calling, which
makes the environment of an `all` run irrelevant; keeping the two invocations separate is still the
only way to get both sets.

Two further regimes exist, are not part of `all` either, and the paper no longer reports them:

```sh
ERROR_NANOPORE=1 sh ./bin/make_fastqs.sh protozoa         # 10 % per-base error, 125 bp
ERROR_NANOPORE_LONG=1 sh ./bin/make_fastqs.sh protozoa    # 10 % per-base error, 3,926 bp
```

InSilicoSeq's `basic` model only substitutes bases, whereas real Nanopore error is indel-heavy, so
these never did more than approximate long-read data -- NanoSim, trained on real Nanopore reads, is
what covers that ground. They are kept for a quick check of how the classification behaves at a high
error rate or a long read length, where training a NanoSim model would be disproportionate, and they
cost hours to generate and to classify, which is why a full run no longer makes them.

Each regime writes its own fastq mapping file, so they never clobber one another. Note that
InSilicoSeq skips every input record shorter than the read length, so the long regime draws only from
the longer contigs and can fall well short of the requested count -- for `protozoa` a request of
2,000 yielded 196. The script reports the shortfall when the result is less than half of what was
asked for.

For `tick-borne`, NanoSim trains an error model on the real Nanopore reads of a tick sample and
applies it to the RefSeq genomes of the twelve tick-borne genera -- one simulated fastq file per
real one, for all eight ticks of the first paper. This needs two things beforehand.

The real reads, from the SRA study PRJNA790938. They are declared by URL in the Genestrip fastq map
`data/fastq/ticks_real.txt`, so Genestrip's own goal `fastqdownload` fetches them -- no download
script of ours is involved:

```sh
mvn exec:exec@fastqdl -Dname=tick-borne -Dfqmap=ticks_real.txt
```

`-ll` puts them in the common fastq directory, `data/fastq`, named after their key in the map:
`tick1.fastq.gz` .. `tick8.fastq.gz`, about 3.2 GB in total. Genestrip skips whatever is already in
place, so the command is safe to repeat. A fastq map has no checksum column, though -- that exists
only on the genome side, in `common/fasta/downloads.txt` -- so nothing verifies what a host actually
sent, which is why `make_fastqs.sh tick-borne` runs `gzip -t` over the files before using them.
`make_fastqs.sh tick-borne` runs it itself for whatever is missing; set `SKIP_FETCH=1` to suppress
that, e.g. when the reads were copied over from another machine.

The map points at NCBI's own SRA fastq endpoint. It serves gzip even though the URL does not say
so, which is what makes the `.fastq.gz` naming correct. If that host is unreachable from your
network, `ticks_real.txt` lists two alternatives (NCBI's other host name, and ENA, which mirrors the
identical runs); switching is a search-and-replace on the URLs, since Genestrip names the download
after the map key rather than after the URL. `make_fastqs.sh tick-borne` verifies every file with
`gzip -t` before starting, so a host that answers with an error page or with uncompressed data fails
immediately instead of feeding NanoSim garbage.

And the genome list

```sh
mvn exec:exec@nanosimlist -Dname=tick-borne
```

which labels every reference genome `<taxid>x<index>`. That label becomes the prefix of every
simulated read's name and is what keeps the read's origin recoverable afterwards. It is derived from
the per-accession FASTA files of the goal `extractrefseqfasta`, so the `tick-borne` database must be
built before it can be created.

Each tick simulation also leaves three files in `results` that describe what it was, none of which
can be recovered afterwards -- NanoSim writes every sample under the same `training_` prefix in one
shared working directory, which is cleared between samples:

| File | What it is |
| --- | --- |
| `tick-borne_<sample>_quantification.tsv` | The abundance table: per reference genome, the share of the simulated reads it supplies. Derived by aligning the *real* reads against the genome list, so it states the real sample's composition as projected onto the twelve genera. |
| `tick-borne_<sample>_error_rate.tsv` | The per-base error NanoSim derived while training, from the alignment rather than from the quality strings -- see below. |
| `tick-borne_nanosim_genomes.tsv` | The genome list itself, which is what resolves an abundance row to a taxon. |

The error rate has to come from the alignment because these runs carry placeholder qualities: the SRA
copies hold two values, Q3 for some 90 % of the bases and Q30 for the rest, NanoSim trains its
base-quality model on exactly those, and an error rate read off the simulated file would therefore be
nonsense. Nothing but the choice of yardstick is affected -- Genestrip classifies on *k*-mers and
ignores quality scores altogether.

Both parts skip whatever is already present and write the fastq mapping files the next step needs.
That skip is worth knowing about for the abundance table: a sample whose `_sim.fastq` is already
there is not simulated again, so no amount of re-running produces the table for it. Delete
`data/fastq/<sample>_sim.fastq` to force the simulation -- which is what step 5 of
`run_all_exps.sh` does, for exactly those samples whose table is missing.

## 4. Real reads, without ground truth

Two collections of real sequencing data are analysed as well, both taken from the first Genestrip
paper: the human saliva runs, matched against `cv`, and the tick samples, matched against `tb`. No
ground truth exists for either, so precision and recall are undefined; what is measured instead is
how far each database variant narrows the species down -- see the paper's section "Estimating the
gain without ground truth".

The tick reads are the same files NanoSim trains on and are already in place after step 3. The
saliva runs are fetched with sra-tools, which `install_tools.sh` installs:

```sh
sh ./bin/fetch_saliva.sh                  # the three runs the first paper used
sh ./bin/fetch_saliva.sh SRR5571991       # a single run
```

**Mind the volume.** These are deep metagenomic runs of 605 to 981 million read pairs, 122 to 198
Gbp each -- roughly 400 GB of gzipped fastq for the three, and transiently about as much again for
the `.sra` files and fasterq-dump's scratch space. Run them one at a time unless the machine has a
spare quarter of a terabyte. Files are named after their accession, exactly as the original
`genestrip-db-exp` project named them, so anything already fetched there is reused rather than
downloaded again.

Then:

```sh
sh ./bin/run_classification_exps.sh real
```

which writes `results/<db>_<key>_specificity.csv`, one row per fastq key:

| Columns | What they hold |
| --- | --- |
| `fastq key`, `sample` | The raw key, and the label the paper prints. Reports join on the raw key, never on the label. |
| `reads`, `classified unrefined`, `classified refined` | Read counts. The two classified counts agree: a refinement changes *where* a read is placed, not whether it is placed. |
| `obs genus only`, `obs genus only share` | The observable subset the measure averages over, and its share of the classified reads. |
| `ungated precision unrefined`, `ungated precision refined` | The measure itself, for each variant. |
| `rho u`, `rho f`, `est prec g u`, `est prec g f` | The calibration factors and the estimated gated precisions that follow from them. |

**This run needs a simulated run first.** Neither precision here can be compared against a gated one,
because no ground truth exists to gate with -- so the two factors are read out of the `rho u` and
`rho f` columns of a simulated run's `_summary.csv`, named by the last argument of `run_real`:
`nanosim` for the ticks, `iss_saliva` for the saliva runs. Run `real` before that summary exists and
the four right-hand columns come out empty, or, worse, are filled from a stale copy still lying in
`results`. A calibration of exactly one row is applied to every sample regardless of its key -- which
is how the single saliva-like read set stands in for all three saliva runs, whose keys are SRA
accessions and match nothing; more than one row means the rows are per sample and only a key match
will do, which is how each tick borrows the calibration derived from its own sample.

The difference between the two ungated precisions is **not** a bound on the precision gain. It
becomes an *estimate* of it only through those factors, which is what the two `est prec g` columns
are: each carries one level from what can be observed to what it stands for, and their difference is
the estimated gain.

### The Kraken tools on the real saliva runs

The real runs have no ground truth, so a comparison with another tool cannot ask whose
classification is right. It can ask something else: of the reads the unrefined `cv` database left at
a genus, where does another tool place those very reads. That is what this adds to Table 9 of the
paper.

```sh
sh ./bin/run_classification_exps.sh kraken     # builds, then classifies simulated AND real
sh ./bin/run_classification_exps.sh real       # evaluates, picking the TSVs up by itself
```

The classification step is `bin/kraken_classify.sh` with a fastq map, which switches it from the
four simulated sets to the runs a map names:

```sh
FQMAP=saliva_real.txt sh ./bin/kraken_classify.sh all
```

Keys and files then come from the map, in its order, so a run that is paired but has only its first
mate listed is read exactly as the Genestrip run reads it.

**The per-read output is filtered as it is produced.** The tool writes into a FIFO and an `awk` reads
it, keeping the three fields the evaluation uses and only the lines whose taxon is not 0. Without
that the output would be 552 GB per tool: the five saliva runs hold 4.1 billion reads and the
simulated files measure 134 bytes per read. Filtered it is a few gigabytes, since only some two per
cent of these reads are classified against a viral database at all. The test is on the taxon rather
than on the `C`/`U` status, which is the same thing in both formats; KrakenUniq would also do it
itself with `--only-classified-output`, but Kraken 2 parses that option and then ignores it, so one
filter serves both.

Dropping the unclassified lines would flatter the tools if the average were taken over what the
evaluation sees, so it is not: `SpecificityReport` divides the sum of scores by the size of the
subset the unrefined run collected, which makes a read the tool never reported count as the zero it
is. The column `obs genus only reported` says how many of the subset each tool did report, so the
gap stays visible.

`SpecificityReport` scores whatever it finds under `results/kraken` for the samples of the run and
writes `<db>_<tool>_<key>_specificity.csv` beside its own. **The subset there is Genestrip's, not the
tool's** — that is the whole point, since both have to be scored on one and the same set of reads.
This is the opposite of the simulated case in `RefinementAccuracyReport`, where each tool collects
its own genus-only subset because the columns of Tables 7 and 8 are per tool. The two differ because
the questions differ, and both classes say so in their Javadoc.

The calibration is the tool's own, from `<db>_<tool>_iss_saliva_summary.csv`. Its `rho f` is empty,
so the refined estimate stays empty as well, as do the refined columns of the row.

## 5. Classification quality experiments

```sh
sh ./bin/run_classification_exps.sh
```

Every simulated fastq file is matched against the unrefined and against the refined database and the
classifications are scored against the reads' known ground truth. The result is one row per fastq
file and database variant in `results/<db>_<key>_accuracy.csv`, holding the boolean positive counts
at the genus and species rank together with the candidate-weighted species count of the paper's
Section "Classification quality" -- the latter being the measure that can distinguish a refined
database from an unrefined one at all, since the boolean counts are blind to a classification that
narrows the species down without reaching a single one.

The column `unresolved` counts reads whose ground truth could not be recovered from their name. It
should be small; a large value means the read identifiers do not fit the accession map, and the
remaining figures then rest on a fraction of the data.

The last seven columns serve the ground-truth-free estimate of the paper's Section "Estimating the
gain without ground truth", and they come in two blocks that must not be crossed:

| Column | What it holds |
| --- | --- |
| `genus only` | The size of the gated subset *R<sub>g</sub>*: the reads the unrefined database placed correctly as far as a genus and no further. Its membership follows from the read's true species, so it cannot be formed on a real sample. |
| `genus only score` | The accumulated candidate weights over that subset. |
| `genus only precision species cand` | Their average, i.e. the restricted precision. Between the two variants its difference is the gain where a gain was possible. |
| `genus only species share` | The share of that subset this variant pins down to a single species -- zero for the unrefined variant by construction. |
| `genus only zero scoring` | How many reads of the subset score zero, i.e. were placed where their true species is not in question at all. These are the whole of the difference between the gated and the ungated measure. |
| `obs genus only` | The size of the *observable* subset *R'<sub>g</sub>*: the reads the unrefined database left at a genus, right or wrong. Its membership follows from the assignment alone. |
| `obs genus only ungated precision` | The one measure of the whole file that never consults the ground truth: the reciprocal number of candidate species averaged over that subset, without testing whether the read's true species is among them. |

The two averages are taken over **different subsets** -- that is the point of having both, and the
easiest thing in the file to get wrong. The gated one divides by `genus only`, the ungated one by
`obs genus only`.

The calibration the paper transfers to a real sample is a **pair** of factors on the levels, not one
factor on the gain: `rho u` carries the unrefined ungated precision to the gated one it stands for
and `rho f` does the same after the refinement. Both are computed in `_summary.csv` (see below),
where the two variants sit in one row and the quotient can be formed at all. They are kept apart
because they answer to different effects, and they are ratios of levels rather than of gains for two
reasons: they stay within [0, 1] and so read as the shrinkage factors they are, where a ratio of two
small differences has no such bound and crosses one on "MiSeq"; and they let a real sample be given
an estimated precision before and after rather than only an estimated gain.

### The files the paper's tables read

`_accuracy.csv` is the raw record, one row per fastq file *and database variant*. That is the right
shape for the measurements and the wrong one for a table, which puts the unrefined and the refined
figure of a measure side by side. Three further files are therefore written beside it, one row per
read set with both variants in it, and it is these the paper includes:

| File | Holds |
| --- | --- |
| `<db>_<key>_summary.csv` | The restricted precision of both variants over both subsets, and the calibration pair `rho u`, `rho f`. The gains get no column: each is the difference of two columns already there. |
| `<db>_<key>_quality.csv` | The boolean counts at the genus and species rank -- given for the *unrefined* variant only, since across every read set none of them moves by more than 0.0007 between the two -- and the candidate-weighted precision, recall and F1 for both. |
| `<db>_<key>_simdata.csv` | What the read sets are rather than how well they classified: read length, per-base error, reads generated, how many fall within the database's scope, and how many were unresolved. |

The first three columns of `_simdata.csv` cannot be recovered at evaluation time -- by then the
settings that produced the reads are gone, and for NanoSim the quality strings are known to be
meaningless. They are measured by `make_fastqs.sh` as it generates each set, left in
`results/<db>_simparams.csv`, and joined in here by fastq key. A read set with no row there still
gets a row, with those columns empty, so that a set generated before the script recorded its
parameters appears as an incomplete row rather than vanishing from the table.

Every one of these files carries the raw `fastq key` in its own column beside the printed label.
Reports join on the raw key, never on the label, so that renaming a label cannot silently break a
join -- and the label is escaped for LaTeX, since a table cell is where it ends up.

### Per-taxon results

The performance part of the run also preserves the match goals' own output, one CSV per fastq key
per database variant, as `results/<goal>_<logkey>_<key>.csv` -- for instance `match_ticks_tick1.csv`
for the real tick 1 against the unrefined database and `ftmatch_tick-borne_tick1.csv` for its NanoSim
simulation against the refined one. Where the reports above aggregate, these say per taxon how many
reads were assigned to it, which is the only per-taxon view of a run that survives it.

The log key in the name is load-bearing: the simulated and the real tick runs use the same keys
`tick1 .. tick8` from two different maps, so in the project folder the second run overwrites the
first's files -- and the `clean` that precedes each measured run deletes them outright.

These ride on the performance measurement, so **no cgmemtime means no per-taxon results**: `run_perf`
returns early without it and the match goals never run.

### The k = 24 control of cv

A second, complete evaluation strand for `cv` alone, on the same genomes, the same taxonomy and the
same four read sets, built at *k* = 24 instead of Genestrip's default of 31:

```sh
sh ./bin/k24_exps.sh all          # the twin, its database, and the four read sets against it
sh ./bin/k24_exps.sh projects     # only the twin project, seconds
sh ./bin/k24_exps.sh build        # only the database and the reports over it
sh ./bin/k24_exps.sh accuracy     # only the classification
```

It is step 11 of `run_all_exps.sh`, so a full run produces it without being asked.

Why *k* = 24. Kraken 2 classifies with a spaced seed of weight 24 — 31 minimizer positions of which
7 are masked — whereas Genestrip matches exact 31-mers. The weight of a seed is the number of
positions that must agree, and under independently placed substitutions it is the whole story: a
spaced seed of weight 24 and a contiguous 24-mer match with the same probability, both against a
read carrying errors and against a related species. A Genestrip database at *k* = 24 therefore sits
where Kraken 2 sits on that trade-off.

It also isolates one factor behind the refinement's gain. *k* is the only knob that moves the share
of *k*-mers stored above the data taxa — the mass a refinement can push down — without touching the
genomes, the read sets, the measures or the code. The strand reports exactly the quantities the main
run reports, under the project name `viral-k24`: the share above the data taxa and the two subtree
precisions from `viral-k24_dbquality.csv` and `viral-k24_ftquality.csv`, and the per-read precisions
from `viral-k24_iss_accuracy.csv`, `viral-k24_iss_perfect_accuracy.csv` and
`viral-k24_iss_saliva_accuracy.csv`.

`bin/k24_projects.sh` writes the twin: symlinks to the original's `taxids.txt`, `additional.txt` and
`categories.txt`, a `config.properties` generated from the original's with `kMerSize=24` appended,
and the fastq maps linked under the twin's own name. The last of these is what makes the two rows
comparable — the twin is scored on the very read files `viral` was scored on, not on a second
simulation of them. `K=27 sh ./bin/k24_exps.sh all` builds `viral-k27` the same way.

Its figures go into CSV files of their own. `CSV_SUFFIX=-k24` keeps the generation timings and disk
sizes in `db_gen_perf-k24.csv` and `db_disk_sizes-k24.csv` rather than overwriting the main run's,
and every other file carries the project name. The twin costs about one to two hours and some 8 GB
of disk beside the original, and it simulates no reads, so nothing here needs InSilicoSeq.

### The sampled control of cv

The second twin of `cv`, and the second half of the same question. `cv_k24` above isolates the seed
weight; this one isolates the store:

```sh
sh ./bin/sampling_exps.sh all          # the twin, its database, and the four read sets against it
sh ./bin/sampling_exps.sh projects     # only the twin project, seconds
```

It runs in step 11 of `run_all_exps.sh` beside the *k* = 24 twin.

Kraken 2 keeps a 4-byte cell per minimizer — 17 bits of truncated MurmurHash3 and 15 bits of taxon
index, nothing of the *k*-mer itself — and enters only the minimizers, which at its defaults *k* = 35
and ℓ = 31 is a window of five and so a density of 2/(5+1) = 1/3. Genestrip enters every *k*-mer in a
64-bit word that holds it exactly, plus a Bloom filter of 10 bits per entry in front. Measured on
`cv` that is 9.25 bytes per entry in memory against 5.71, a factor no sampling rate can change, so a
Genestrip database of Kraken 2's size needs a rate of about four. `kMerSampling=4` is what the twin
`viral-s4` sets.

What it costs is paid on the reads. The sampling selects by the *k*-mer and not by its position
(`KMerSampling.java` uses a multiplicative threshold, `kmer * 0x9E3779B97F4A7C15` unsigned against
`(2^64-1)/n`), which is what keeps the tax ids right: a *k*-mer is entered in every genome it occurs
in or in none of them. A read therefore keeps about one *k*-mer in four at random positions. On clean
reads that is harmless; on error-rich ones the surviving error-free *k*-mers are thinned by the same
factor, and that is the number this strand produces.

`CSV_SUFFIX=-s4` keeps its timings and disk sizes in `db_gen_perf-s4.csv` and
`db_disk_sizes-s4.csv`. `S=5 sh ./bin/sampling_exps.sh all` builds `viral-s5` instead.

`K` sets a *k*-mer length along with the sampling, so

```sh
K=24 sh ./bin/sampling_exps.sh all     # viral-k24-s4
```

builds the twin that is smaller both ways at once. The two halvings act on different terms: the seed
weight on how many *k*-mers are distinct, the sampling on how many of them are entered. Only the
combination reaches Kraken 2's size — on `cv` 671 MiB against roughly 702 — and it is the only one of
the three whose seed weight and whose entry count are both comparable with Kraken 2's. The twin is
always built from the original project, never from the `-k24` twin, so one project folder holds one
complete configuration.

## 6. Machine description

The paper states what hardware the experiments ran on, mirroring the corresponding paragraph of the
first Genestrip paper. Run this **on the machine that executed the experiments**, ideally right
after a full run:

```sh
sudo sh ./bin/sysinfo.sh
```

It writes two files:

| File | Purpose |
| --- | --- |
| `results/sysinfo.txt` | The full report -- CPU, memory, DMI/BIOS, storage, PCI, OS, every installed JVM. For the record and for answering reviewer questions. |
| `results/sysinfo.tex` | Just the facts the paper quotes, as LaTeX macros. |

`sudo` only affects the report: DMI, memory-slot and SMART details are root-only. Everything in the
LaTeX file comes from `lscpu`, `/proc`, `lsblk` and `java`, so it is complete without it.

`paper.tex` inputs `results/sysinfo.tex` and uses `\sysRamGB`, `\sysCores`, `\sysWorkers`,
`\sysDiskSize`, `\sysDiskType`, `\sysCpuModel`, `\sysOs`, `\sysJavaVendor`, `\sysJavaVersion` and
`\sysJavaVm` in Section "Database refinement performance".

`\sysWorkers` is the number of worker threads the goals actually run with. `pom.xml` sets
`gs.threads` to `-1`, which Genestrip reads as one thread per available processor less one, so the
macro is simply `\sysThreads - 1`. That is also what the classification experiments use, since
`AccuracyEvaluator` sets the same value in code and the `accuracy` execution passes no thread
option -- keeping the property at `-1` is what makes every experiment of the paper run with the
same thread count. Pin it with `-Dgs.threads=19` to reproduce a specific run; the paper's figure
then has to be set by hand, since the macro follows the hardware rather than the override.

Two safeguards make a missing or incomplete file obvious instead of silent:

* Anything the script cannot determine becomes `\sysUnknown`, which typesets as a bold **??**.
* If the file is absent altogether, `paper.tex` falls back to the same marker, so it still compiles.

So if the PDF shows **??** anywhere in that paragraph, the script has not been run on the right
machine yet, or a value needs to be filled in by hand. The script also lists what it could not
determine, both on the terminal and at the end of `sysinfo.txt`.

Note that the JVM reported is the one on `PATH` at the time the script runs, which is only the one
that executed the goals if the environment is the same -- so run it from the same shell as the
experiments.

## Results

All experiments write to `./results`. Those files are meant to be consumed as they are: a paper
includes the CSV, SVG and LaTeX fragments by name from a results folder of its own, so copy them
there. `run_all_exps.sh` does that copy as its last step when `PAPER_RESULTS` names the folder, and
otherwise leaves everything in `./results`.

**Copying the folder by hand means one more command.** Two of the files a paper reads are derived
from the others rather than written by a goal: `dbstats.tex`, which holds every figure the text
states as a macro, and `matchperf.csv`, which is condensed from the `perf_*.log` files. The copy
step regenerates the first in its target; an rsync does not. So after copying, run

```sh
sh ./bin/paper_stats.sh /path/to/the/paper/results
```

against the folder you copied into — otherwise its tables keep the numbers of the previous batch
while the CSVs beside them are new, which is exactly the mismatch that script exists to prevent.
`matchperf.csv` needs nothing extra: step 12 writes it into `./results` before the copy.

| File | Written by | Section |
| --- | --- | --- |
| `<db>_dbinfo.csv`, `<db>_ftdbinfo.csv` | `run_exps.sh` | 2 |
| `<db>_dbquality.csv`, `<db>_ftquality.csv` | `run_exps.sh` | 2 |
| `<db>_kmerrankstatscsv.csv`, `<db>_branchhistorankcsv.csv` | `run_exps.sh` | 2 |
| `db_gen_<db>.log`, `ftdb_gen_<db>.log` | `run_exps.sh` | 2 |
| `<db>_simparams.csv` | `make_fastqs.sh` | 3 |
| `tick-borne_<sample>_quantification.tsv`, `_error_rate.tsv`, `tick-borne_nanosim_genomes.tsv` | `make_fastqs.sh` | 3 |
| `<db>_<key>_accuracy.csv`, `_summary.csv`, `_quality.csv`, `_simdata.csv` | `run_classification_exps.sh` | 5 |
| `<db>_<key>_specificity.csv` | `run_classification_exps.sh` | 4 |
| `match_<logkey>.log`, `ftmatch_<logkey>.log` | `run_classification_exps.sh` | 5 |
| `<goal>_<logkey>_<key>.csv` | `run_classification_exps.sh` | 5 |
| `perf_<scenario>_<key>_<goal>.log`, `matchperf.csv` | `perf_scenarios.sh` | 5 |
| `viral-k24_*`, `db_gen_perf-k24.csv`, `db_disk_sizes-k24.csv` | `k24_exps.sh` | 5 |
| `viral-s4_*`, `db_gen_perf-s4.csv`, `db_disk_sizes-s4.csv` | `sampling_exps.sh` | 5 |
| `viral-k24-s4_*`, `db_gen_perf-k24-s4.csv`, `db_disk_sizes-k24-s4.csv` | `K=24 sampling_exps.sh` | 5 |
| `dbstats.tex` | `paper_stats.sh` | below |
| `sysinfo.txt`, `sysinfo.tex` | `sysinfo.sh` | 6 |

Among them is `dbstats.tex`, written by `./bin/paper_stats.sh` as the last step of `run_exps.sh`. It
holds the numbers of the paper's database tables -- stored *k*-mers, taxa at the species rank and
below, and the subtree precision per genus -- as LaTeX macros, so that they are read from the CSVs
instead of being copied by hand. Run the script on its own after regenerating individual CSVs:

```sh
sh ./bin/paper_stats.sh
```

Copy it over **together with** the CSVs it was derived from. A `dbstats.tex` beside a different set
of CSVs is the one inconsistency it cannot detect itself; a missing or incomplete one is harmless,
since the paper renders every value it cannot find as a bold `??`.

## Cleaning up

```sh
sh ./bin/clean_all.sh
```

removes the generated databases of every project. The downloaded RefSeq and Genbank data survives,
since the corresponding goals are excluded from the recursive clean -- re-downloading it takes far
longer than rebuilding the databases from it.
