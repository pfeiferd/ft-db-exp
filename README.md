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

Several projects list the human genome in their `additional.txt`, `viral` among them since `cv` is
applied to human saliva. The taxid there tells the LCA update whose genome the file is; it must
**not** also appear in `taxids.txt`, or the database would store the human genome rather than merely
account for it. Without that entry every *k*-mer a virus shares with the human genome -- endogenous
retroviral sequence, integrated herpesvirus, host contamination in the viral assemblies -- stays
stored under its virus and turns human reads into viral hits, which on a saliva sample is most of
the data. The file itself is declared once in `data/common/fasta/downloads.txt` and shared.

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

`all` generates the two Illumina sets and the error-free one of every InSilicoSeq project, which is
what the paper reports. "MiSeq" and "HiSeq" differ in read length *and* error profile at once, so the
error-free set is what separates the two effects: it keeps the read length and drops the errors.

```sh
sh ./bin/make_fastqs.sh protozoa                          # MiSeq (301 bp) and HiSeq (126 bp)
ERROR_FREE=1 sh ./bin/make_fastqs.sh protozoa             # "perfect": same length, no errors
```

Two further regimes exist but are **not** part of `all`, and the paper no longer reports them:

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

Both parts skip whatever is already present and write the fastq mapping files the next step needs.

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

which writes `results/<db>_<key>_specificity.csv`: one row per fastq key with the size of the
observable subset, the ungated precision of each database variant and their difference.

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

The last six columns serve the ground-truth-free estimate of the paper's Section "Estimating the
gain without ground truth". `obs genus only ungated precision` is the one measure of the whole file
that never consults the ground truth: it averages the reciprocal number of candidate species over
the reads the unrefined database left at a genus, without testing whether the read's true species is
among them. Its delta between the two variants is the *specificity* gain, an upper bound on the
precision gain, and the quotient of the two -- computable only here, where the truth is known -- is
the calibration factor the paper calls rho_d. `obs genus only also true` divided by `obs genus only`
says how faithfully the observable subset reproduces the real one; it is essentially 1 on simulated
reads and cannot be expected to be on a real sample. `genus only gate missed` counts the reads that
make rho_d fall short of one.

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
