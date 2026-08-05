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

Individual steps can be run on their own, e.g. only the two databases of one project:

```sh
mvn exec:exec@db -Dname=viral -Dgoal=db      # unrefined
mvn exec:exec@db -Dname=viral -Dgoal=ftdb    # refined, on top of the unrefined one
```

## 3. Simulated reads

```sh
sh ./bin/make_fastqs.sh              # both parts
sh ./bin/make_fastqs.sh viral        # InSilicoSeq only
sh ./bin/make_fastqs.sh tick-borne   # NanoSim only
```

For `viral`, InSilicoSeq applies its Illumina "MiSeq" and "HiSeq" error models to all RefSeq genomes
of the category "Viral", one million reads each as in the first paper. `N_READS=10k` produces a much
smaller set for a quick smoke test.

For `tick-borne`, NanoSim trains an error model on the real Nanopore reads of a tick sample and
applies it to the RefSeq genomes of the twelve tick-borne genera. This needs two things beforehand:
the real reads in `data/fastq` (the script's header lists the SRA accessions) and the genome list

```sh
mvn exec:exec@nanosimlist -Dname=tick-borne
```

which labels every reference genome `<taxid>x<index>`. That label becomes the prefix of every
simulated read's name and is what keeps the read's origin recoverable afterwards.

Both parts skip whatever is already present and write the fastq mapping files the next step needs.

## 4. Classification quality experiments

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

## 5. Machine description

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

All experiments write to `./results`. Those files are consumed directly by the paper: copy them into
its own `results` folder, i.e. `genestrip-docs2/ft-paper/results`, from where the LaTeX sources
include the CSV, SVG and LaTeX fragments by name.

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
