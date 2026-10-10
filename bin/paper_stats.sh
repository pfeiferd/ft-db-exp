#!/bin/sh
#
# Extracts the database statistics the paper's tables state into LaTeX macros, so that they are read
# from the result CSVs instead of being copied by hand. Writes <results>/dbstats.tex.
#
# Usage:
#   sh ./bin/paper_stats.sh [RESULTS_DIR]
#
# RESULTS_DIR defaults to ./results. Copy the generated file into the paper's own results folder
# together with the CSVs it was derived from -- the two belong together, and a stale dbstats.tex
# next to fresh CSVs is exactly the mismatch this script exists to prevent.
#
# Which columns are used, and why it matters:
#
#   stored k-mers      the `stored kmers' field of the TOTAL row of <db>_dbinfo.csv. A refinement
#                      reassigns k-mers but adds none, so the same number describes both variants.
#
#   species and below  the number of rows in <db>_dbinfo.csv whose rank is `species' or one of the
#                      ranks beneath it. The synthetic ranks Genestrip uses to mark where a k-mer
#                      came from -- DATA, REFINED, FILE, ID -- are not taxonomic and are excluded,
#                      as are rows with no rank of their own.
#
#   subtree precision  the `subtree precision' column of <db>_dbquality.csv and <db>_ftquality.csv.
#                      NOT `weighted avg precision', which sits right beside it and is a different
#                      quantity: that one is the quotient of the summed counts, this one averages
#                      the per-k-mer precision p(a) over the k-mers of the subtree. Only the latter
#                      is the sp(m) that the paper defines. The difference is not cosmetic -- for
#                      the viral root it is 0.05 against 0.92 -- because a ratio of sums is
#                      dominated by the few k-mers sitting at nodes with hundreds of data taxa
#                      underneath them.
#
#   restricted sp      the same average taken over the k-mers stored *above* the data taxa only.
#                      After the leaf transformation a k-mer at a data taxon has |D| = 1 and hence
#                      p(a) = 1 by construction, and such k-mers are the bulk of a database -- 99 %
#                      of parasites, 62 % of tick-borne -- so they dominate sp and can never be
#                      improved. Genestrip reports it as `restricted subtree precision', computed in
#                      the same pass as sp itself, so numerator and denominator agree by
#                      construction. It is taken from there and never recomputed here.
#
set -e

scriptdir=$(dirname "$0")
cd "$scriptdir/.."
results=${1:-./results}

[ -d "$results" ] || { echo "No such directory: $results" >&2; exit 1; }

python3 - "$results" <<'PY'
import csv, os, sys

results = sys.argv[1]
out = os.path.join(results, 'dbstats.tex')

# The ranks at or below the species rank. Everything else is either higher up or one of Genestrip's
# synthetic markers (DATA, REFINED, FILE, ID), which say where a k-mer came from rather than what
# the organism is.
SPECIES_AND_BELOW = {'species', 'subspecies', 'strain', 'isolate', 'serotype', 'serogroup',
                     'genotype', 'biotype', 'morph', 'forma', 'forma specialis', 'varietas',
                     'subvariety', 'pathogroup'}

# The last two are the twins of `viral': `viral-k24' from bin/k24_exps.sh is built at k = 24, and
# `viral-k24-s4' from `K=24 bin/sampling_exps.sh' enters one 24-mer in four besides, which is the
# configuration that comes out the size Kraken 2's database is. They are listed here so that their
# share above the data taxa and their two subtree precisions can be read from the same macros as
# every other database's. A twin that was never built has no dbinfo.csv and is skipped, as any
# absent database is.
DBS = ['viral', 'tick-borne', 'strepto', 'nocardia', 'protozoa', 'parasites', 'orthopox',
       'viral-k24', 'viral-k24-s4']

def subtree_data_kmers(info):
    """For every node of a dbinfo file, the k-mers stored at data taxa in its subtree.

    The file lists the tree in preorder with an explicit level, so a node's subtree runs until the
    next row whose level is not greater than its own. Returned per tax id, which is what the
    quality files key on.
    """
    nodes = []
    for r in info:
        try:
            level = int((r.get('level') or '').strip())
            kmers = int((r.get('stored kmers') or '0').strip())
        except ValueError:
            continue
        nodes.append((level, (r.get('taxid') or '').strip(),
                      kmers if (r.get('rank') or '').strip() == 'DATA' else 0))
    out = {}
    for i, (level, taxid, _) in enumerate(nodes):
        total = 0
        for j in range(i, len(nodes)):
            if j > i and nodes[j][0] <= level:
                break
            total += nodes[j][2]
        if taxid:
            out[taxid] = total
    return out

def rows(path):
    if not os.path.exists(path):
        return []
    with open(path, encoding='utf-8') as fh:
        return list(csv.DictReader(fh, delimiter=';'))

def emit(fh, kind, key, value):
    fh.write('\\expandafter\\def\\csname dbstat@%s@%s\\endcsname{%s}\n' % (kind, key, value))

entries = missing = 0
collisions = []
withoutcolumn = set()
aboveshare = {}
with open(out, 'w', encoding='utf-8') as fh:
    fh.write('% Database statistics for the paper\'s tables, generated by ft-db-exp2/bin/paper_stats.sh\n'
             '% from the CSV files beside this one. Do not edit: rerun the script instead.\n'
             '% Look up a value with \\dbstat{<kind>}{<key>}, e.g. \\dbstat{kmers}{viral} or\n'
             '% \\dbstat{spf}{protozoa/Plasmodium}. An unknown key renders as a visible marker.\n')
    for db in DBS:
        info = rows(os.path.join(results, db + '_dbinfo.csv'))
        if info:
            total = [r for r in info if (r.get('name') or '').strip() == 'TOTAL']
            if total:
                emit(fh, 'kmers', db, (total[0].get('stored kmers') or '').strip())
                entries += 1
            n = sum(1 for r in info if (r.get('rank') or '').strip() in SPECIES_AND_BELOW)
            emit(fh, 'species', db, str(n))
            entries += 1
            # The species rank alone, without the strains and subspecies the count above includes.
            # This is the size of the candidate set a read left at the genus of a one-genus database
            # still has open, so it is the reciprocal of that database's ungated precision there.
            n = sum(1 for r in info if (r.get('rank') or '').strip() == 'species')
            emit(fh, 'speciesrank', db, str(n))
            entries += 1
        else:
            missing += 1
        # Only the genera and the root are emitted, and they are keyed by name. Taxon names are not
        # unique across ranks -- protozoa has both a genus and a subgenus called Plasmodium, with
        # different values -- so anything else would silently overwrite. A collision within the
        # emitted set is reported rather than resolved by whichever row happens to come last.
        for kind, suffix, infosuffix in (('spu', '_dbquality.csv', '_dbinfo.csv'),
                                         ('spf', '_ftquality.csv', '_ftdbinfo.csv')):
            leaves = subtree_data_kmers(rows(os.path.join(results, db + infosuffix)))
            seen = {}
            for r in rows(os.path.join(results, db + suffix)):
                v = (r.get('subtree precision') or '').strip()
                name = (r.get('name') or '').strip()
                rank = (r.get('rank') or '').strip()
                if rank != 'genus' and name != 'root':
                    continue
                if not v or v.lower() in ('nan', 'null') or not name.replace(' ', '').isalnum():
                    continue
                key = '%s/%s' % (db, name.replace(' ', ''))
                if key in seen and seen[key] != v:
                    print('  COLLISION %s: %s and %s - table would show whichever came last'
                          % (key, seen[key], v))
                    collisions.append(key)
                    continue
                if key in seen:
                    continue
                seen[key] = v
                emit(fh, kind, key, v)
                entries += 1
                # The same average over the k-mers above the data taxa, where a refinement can act.
                taxid = (r.get('taxid') or '').strip()
                try:
                    K = int((r.get('subtree kmers') or '0').strip())
                except ValueError:
                    K = 0
                L = leaves.get(taxid)
                above = (r.get('subtree kmers above data') or '').strip()
                # The share of the subtree's k-mers that sit above the data taxa, i.e. the part a
                # refinement can still move and the one restricted sp averages over. A refinement
                # moves k-mers onto nodes it enters itself, never onto a data taxon, so the share
                # is the same in both variants; the second pass only checks that.
                if above.isdigit() and K > 0:
                    share = '%.4f' % (100.0 * int(above) / K)
                    if key not in aboveshare:
                        aboveshare[key] = share
                        emit(fh, 'above', key, share)
                        entries += 1
                    elif aboveshare[key] != share:
                        print('  COLLISION %s: k-mers above the data taxa differ between the '
                              'variants, %s against %s' % (key, aboveshare[key], share))
                        collisions.append(key)
                restricted = (r.get('restricted subtree precision') or '').strip()
                if restricted and restricted.lower() not in ('nan', 'null'):
                    emit(fh, kind + 'x', key, restricted)
                    entries += 1
                    if name == 'root' and above.isdigit() and K > 0:
                        emit(fh, 'leafshare', db, '%.3f' % (100.0 * (K - int(above)) / K))
                        entries += 1
                elif above == '0':
                    # Undefined rather than missing: every k-mer of this subtree sits at a data
                    # taxon, where p(a) = 1 by construction, so sp* has nothing to average over.
                    # An em-dash says so; the marker below would claim the measurement failed.
                    # Braced, because the column is an siunitx S column and this is text.
                    emit(fh, kind + 'x', key, '{---}')
                    entries += 1
                else:
                    # No fallback. Deriving the value from the dbinfo file instead means mixing two
                    # accountings: the quality file counts only the k-mers for which p is defined,
                    # dbinfo counts every stored k-mer. At genus level the data-taxon share can then
                    # exceed the subtree total, the subtraction turns negative, and values such as
                    # -2902 result. A missing column is reported and shows as a marker in the paper;
                    # a wrong number would not be noticed.
                    withoutcolumn.add(db)

    # The real-read runs. One <db>_<key>_specificity.csv per read collection, one row per sample,
    # and for `strepto' eighty-three of them -- too many to print, and a mean over them is a
    # computation rather than a lookup, so it belongs here and not in the LaTeX.
    #
    # The average is weighted by |R'_g|, the observable subset each row is itself an average over.
    # Weighting by sample instead would let a sample with one such read count as much as one with
    # four thousand, and `strepto' has both.
    def specrows(path):
        out = []
        for r in rows(path):
            try:
                n = int((r.get('obs genus only') or '0').strip())
                reads = int((r.get('reads') or '0').strip())
                cls = int((r.get('classified unrefined') or '0').strip())
                pu = (r.get('ungated precision unrefined') or '').strip()
                pf = (r.get('ungated precision refined') or '').strip()
                # Present only where a calibration was supplied; SpecificityReport leaves them
                # empty otherwise, and a group then reports no estimate rather than a wrong one.
                eu = (r.get('est prec g u') or '').strip()
                ef = (r.get('est prec g f') or '').strip()
            except ValueError:
                continue
            if n <= 0 or not pu or not pf:
                continue
            out.append((r.get('sample', '').strip(), n, float(pu), float(pf), reads, cls,
                        float(eu) if eu else None, float(ef) if ef else None))
        return out

    # Everything a row of Table \ref{realgain} needs, for one group of samples. `strepto' prints
    # groups rather than samples: eighty-three rows would not fit, and the split by what the
    # reference standard says about the pneumococcus is the comparison the case study makes.
    def emitgroup(fh, db, tag, sel):
        if not sel:
            return 0
        N = sum(x[1] for x in sel)
        pu = sum(x[1] * x[2] for x in sel) / N
        pf = sum(x[1] * x[3] for x in sel) / N
        reads = sum(x[4] for x in sel)
        cls = sum(x[5] for x in sel)
        emit(fh, 'rg', '%s/%ssamples' % (db, tag), str(len(sel)))
        emit(fh, 'rg', '%s/%ssubset' % (db, tag), str(N))
        # The smallest and the largest observable subset among the samples of the group. The prose
        # states that range, and states it from here rather than from a reading of the table, so
        # that a collection gaining a sample does not leave a stale pair of numbers behind.
        emit(fh, 'rg', '%s/%ssubsetmin' % (db, tag), str(min(x[1] for x in sel)))
        emit(fh, 'rg', '%s/%ssubsetmax' % (db, tag), str(max(x[1] for x in sel)))
        emit(fh, 'rg', '%s/%sgroupreads' % (db, tag), str(reads))
        emit(fh, 'rg', '%s/%sgroupclassified' % (db, tag), str(cls))
        emit(fh, 'rg', '%s/%sgroupshare' % (db, tag), '%.2f' % (100.0 * N / cls) if cls else '0')
        emit(fh, 'rg', '%s/%sprecu' % (db, tag), '%.4f' % pu)
        emit(fh, 'rg', '%s/%sprecf' % (db, tag), '%.4f' % pf)
        emit(fh, 'rg', '%s/%sgain' % (db, tag), '%.0f' % (100.0 * (pf - pu) / pu) if pu else '0')
        emit(fh, 'rg', '%s/%sopenu' % (db, tag), '%.1f' % (1.0 / pu) if pu else '0')
        emit(fh, 'rg', '%s/%sopenf' % (db, tag), '%.1f' % (1.0 / pf) if pf else '0')
        # The gated estimates, averaged over |R'_g| exactly as the precisions above are. Emitted
        # only when every sample of the group carries one, since an average over part of a group
        # would be read as an average over all of it.
        n_emitted = 10
        if all(x[6] is not None and x[7] is not None for x in sel):
            emit(fh, 'rg', '%s/%sestu' % (db, tag), '%.4f' % (sum(x[1] * x[6] for x in sel) / N))
            emit(fh, 'rg', '%s/%sestf' % (db, tag), '%.4f' % (sum(x[1] * x[7] for x in sel) / N))
            n_emitted += 2
        return n_emitted

    # (database, report key, macro tag). The tag keeps two collections of one database apart:
    # `strepto' is run on its clinical samples and on the saliva runs, and both write their figures
    # under the same database name. An empty tag is the collection the prose speaks of by default.

    for db, key, tag in (('strepto', 'saliva', 'saliva'), ('tick-borne', 'ticks', ''),
                         ('viral', 'saliva', '')):
        path = os.path.join(results, '%s_%s_specificity.csv' % (db, key))
        allrows = rows(path)
        if not allrows:
            continue
        sel = specrows(path)
        totreads = sum(int((r.get('reads') or '0').strip()) for r in allrows)
        totcls = sum(int((r.get('classified unrefined') or '0').strip()) for r in allrows)
        emit(fh, 'rg', '%s/%sreads' % (db, tag), str(totreads))
        emit(fh, 'rg', '%s/%sclassified' % (db, tag), str(totcls))
        emit(fh, 'rg', '%s/%sallsamples' % (db, tag), str(len(allrows)))
        entries += 3
        n = emitgroup(fh, db, tag, sel)
        entries += n
        if sel and totcls:
            share = 100.0 * sum(x[1] for x in sel) / totcls
            emit(fh, 'rg', '%s/%sshare' % (db, tag), '%.2f' % share)
            emit(fh, 'rg', '%s/%simproved' % (db, tag), str(sum(1 for x in sel if x[3] > x[2])))
            emit(fh, 'rg', '%s/%sflat' % (db, tag), str(sum(1 for x in sel if x[3] == x[2])))
            emit(fh, 'rg', '%s/%sworse' % (db, tag), str(sum(1 for x in sel if x[3] < x[2])))
            entries += 4

print('wrote %s with %d entries' % (out, entries))
if missing:
    print('  %d database(s) had no dbinfo.csv in %s and were skipped' % (missing, results))
if collisions:
    print('  %d collision(s) - those keys were left out entirely' % len(collisions))
if withoutcolumn:
    print("  no 'restricted subtree precision' column in: %s" % ', '.join(sorted(withoutcolumn)))
    print('  -> regenerate those quality CSVs with a current Genestrip-FT; the paper shows a marker')
PY
