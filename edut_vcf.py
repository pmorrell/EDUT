#!/usr/bin/env python3
"""
edut_vcf.py - VCF front end for the EDUT "pattern a" triplet error-detection test.

Ports the core algorithm from EDUT_1.2.1.pl (Toleno et al. 2009) to Python 3
and adapts the input stage to read biallelic SNPs directly from a VCF instead
of an aligned FASTA or ms-style simulation file.

Core idea, unchanged from the Perl version:
  - Encode each haplotype at each segregating site as 0 (matches the more
    common/ancestral-like state -- here, REF) or 1 (ALT).
  - For every ordered triplet of sites A < B < C, run the four-gamete test on
    the pairs (A,C) and (A,B) and (B,C). "Pattern a" is: A/C compatible
    (only 3 of 4 gametes observed) while A/B and B/C are each incompatible
    (all 4 gametes observed). This pattern is unlikely to arise from a single
    recombination event and instead points at a genotyping or phasing error
    at site B (or at one of the rare haplotypes involved).
  - The haplotype(s) carrying the rare (count <= 1) gamete class in each
    incompatible pair are flagged and accumulate a raw priority score per
    site. The raw score is corrected by the overall frequency of pattern a
    and by a "middle site" correction (a site in the middle of the window
    has more chances to be flagged than one at either end).

Phasing:
  The four-gamete/triplet test is only meaningful on haplotypes (gametic
  phase known). If every genotype in the window is homozygous (0|0, 1|1, or
  equivalently 0/0 or 1/1 -- the GT separator doesn't matter for a
  homozygous call), each diploid individual is collapsed to a single
  haploid sequence. Phasing is therefore NOT required for inbred/homozygous
  data, whether the VCF itself is phased or not, and flags still trace back
  to one individual. As soon as a heterozygous genotype appears, unambiguous
  phase (the '|' GT separator) is required to build haplotypes; an unphased
  heterozygous genotype is reported by SNP position and sample name and
  treated as missing at that site, with a warning that this may result in
  unidentified genotyping errors.
"""

import argparse
import gzip
import itertools
import math
import sys
from collections import defaultdict

CUTOFF_FREQ = 1  # a gamete-class count at or below this is "rare" and gets flagged


def open_maybe_gzip(path):
    if path == "-":
        return sys.stdin
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


class Site:
    __slots__ = ("chrom", "pos", "ref", "alt", "alleles")

    def __init__(self, chrom, pos, ref, alt, alleles):
        self.chrom = chrom
        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.alleles = alleles  # per-haplotype '0'/'1'/'N', length = n_haplotypes


def parse_gt(gt_field):
    """Return (allele1, allele2, phased) as strings ('0','1','.') from a GT value."""
    sep = "|" if "|" in gt_field else "/"
    phased = sep == "|"
    parts = gt_field.split(sep)
    if len(parts) == 1:
        # haploid call (e.g. male sex chromosome, or already-haploid VCF)
        return parts[0], None, phased
    return parts[0], parts[1], phased


def read_vcf_sites(path, max_missing, exclude_near_indel, allow_unphased):
    """Yield Site objects for biallelic SNPs, plus sample names.

    Each haplotype is tracked separately (2 per diploid sample, in VCF
    column order) so downstream windowing can decide how to collapse them.
    """
    samples = []
    indel_positions_by_chrom = defaultdict(list)
    raw_records = (
        []
    )  # keep everything so we can filter by proximity to indels after one pass

    with open_maybe_gzip(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                fields = line.split("\t")
                samples = fields[9:]
                continue
            fields = line.split("\t")
            chrom, pos, _id, ref, alt_field = (
                fields[0],
                int(fields[1]),
                fields[2],
                fields[3],
                fields[4],
            )
            fmt = fields[8].split(":")
            gt_index = fmt.index("GT") if "GT" in fmt else None
            alts = alt_field.split(",")
            is_indel = len(ref) != 1 or any(
                len(a) != 1 for a in alts if a not in (".", "*")
            )
            if is_indel:
                indel_positions_by_chrom[chrom].append(pos)
                continue
            if len(alts) != 1 or alts[0] in (".", "*"):
                # skip multiallelic / no-call sites; EDUT expects biallelic 0/1 encoding
                continue
            if gt_index is None:
                continue
            raw_records.append((chrom, pos, ref, alts[0], fields[9:], gt_index))

    for chrom, pos, ref, alt, sample_fields, gt_index in raw_records:
        if exclude_near_indel > 0:
            near = any(
                abs(pos - ip) <= exclude_near_indel
                for ip in indel_positions_by_chrom.get(chrom, [])
            )
            if near:
                continue

        haplotype_alleles = []
        n_missing_haplotypes = 0
        for sample_idx, sample_field in enumerate(sample_fields):
            subfields = sample_field.split(":")
            gt = subfields[gt_index] if gt_index < len(subfields) else "./."
            a1, a2, phased = parse_gt(gt)
            # phase is only meaningful for heterozygous calls; an unphased
            # het can't be assigned to a haplotype without guessing
            if a2 is not None and a1 != a2 and not phased and not allow_unphased:
                sample_name = (
                    samples[sample_idx]
                    if sample_idx < len(samples)
                    else f"sample{sample_idx}"
                )
                print(
                    f"# WARNING: unphased heterozygous genotype at {chrom}:{pos} in sample "
                    f"'{sample_name}' treated as missing; this may result in unidentified "
                    f"genotyping errors.",
                    file=sys.stderr,
                )
                haplotype_alleles.append("N")
                haplotype_alleles.append("N")
                n_missing_haplotypes += 2
                continue
            for allele in (a1, a2):
                if allele is None:
                    continue
                if allele == ".":
                    haplotype_alleles.append("N")
                    n_missing_haplotypes += 1
                elif allele == "0":
                    haplotype_alleles.append("0")
                elif allele == "1":
                    haplotype_alleles.append("1")
                else:
                    # third+ allele observed in this sample though ALT is biallelic overall
                    haplotype_alleles.append("N")
                    n_missing_haplotypes += 1

        if not haplotype_alleles:
            continue
        missing_fraction = n_missing_haplotypes / len(haplotype_alleles)
        if missing_fraction > max_missing:
            continue
        # drop monomorphic/singleton sites: not informative for the four-gamete test
        n_derived = sum(1 for a in haplotype_alleles if a == "1")
        n_called = sum(1 for a in haplotype_alleles if a != "N")
        if n_derived < 2 or (n_called - n_derived) < 2:
            continue

        yield Site(chrom, pos, ref, alt, haplotype_alleles)

    return samples


def build_haplotype_labels(samples):
    labels = []
    for s in samples:
        labels.append(f"{s}_hap1")
        labels.append(f"{s}_hap2")
    return labels


def collapse_to_haploid_if_possible(sites, samples):
    """If every genotype in the window is homozygous, collapse each sample's
    two haplotypes into one haploid call. Returns (sites_maybe_collapsed,
    labels, collapsed: bool). Raises ValueError if heterozygous + unphased
    genotypes are mixed such that phase is required but unavailable --
    caller should have already turned unresolvable hets into 'N'.
    """
    all_homozygous = True
    for site in sites:
        for i in range(0, len(site.alleles), 2):
            a1, a2 = site.alleles[i], site.alleles[i + 1]
            if a1 != "N" and a2 != "N" and a1 != a2:
                all_homozygous = False
                break
        if not all_homozygous:
            break

    if not all_homozygous:
        return sites, build_haplotype_labels(samples), False

    collapsed_sites = []
    for site in sites:
        collapsed_alleles = []
        for i in range(0, len(site.alleles), 2):
            a1, a2 = site.alleles[i], site.alleles[i + 1]
            if a1 == "N" and a2 == "N":
                collapsed_alleles.append("N")
            elif a1 == "N":
                collapsed_alleles.append(a2)
            elif a2 == "N":
                collapsed_alleles.append(a1)
            else:
                collapsed_alleles.append(a1)  # a1 == a2, homozygous
        collapsed_sites.append(
            Site(site.chrom, site.pos, site.ref, site.alt, collapsed_alleles)
        )
    return collapsed_sites, list(samples), True


def four_gamete(col_a, col_b):
    """Return (n00, n10, n01, n11) if all four gametic classes are present
    (ignoring haplotypes missing at either site), else None."""
    n00 = n10 = n01 = n11 = 0
    for a, b in zip(col_a, col_b):
        if a == "N" or b == "N":
            continue
        if a == "0" and b == "0":
            n00 += 1
        elif a == "1" and b == "0":
            n10 += 1
        elif a == "0" and b == "1":
            n01 += 1
        else:
            n11 += 1
    if n00 and n10 and n01 and n11:
        return (n00, n10, n01, n11)
    return None


def find_pattern_a(sites, max_span):
    """Return list of dicts describing each 'pattern a' triplet found, plus
    the total triplet count considered (for the proportion statistic)."""
    n_sites = len(sites)
    columns = [site.alleles for site in sites]

    # precompute pairwise four-gamete results
    lookup = {}
    for i, j in itertools.combinations(range(n_sites), 2):
        lookup[(i, j)] = four_gamete(columns[i], columns[j])

    results = []
    n_triplets_considered = 0
    for a, b, c in itertools.combinations(range(n_sites), 3):
        if max_span is not None and (sites[c].pos - sites[a].pos) > max_span:
            continue
        n_triplets_considered += 1
        ac = lookup[(a, c)]
        ab = lookup[(a, b)]
        bc = lookup[(b, c)]
        if ac is None and ab is not None and bc is not None:
            results.append({"a": a, "b": b, "c": c, "ab_counts": ab, "bc_counts": bc})
    return results, n_triplets_considered


def flag_haplotype(rare_pair, pos1, pos2, columns):
    """Find the haplotype index carrying the rare 2-character gamete class
    ('00', '10', '01', or '11') at the two given site indices."""
    for hap_idx in range(len(columns[pos1])):
        test = columns[pos1][hap_idx] + columns[pos2][hap_idx]
        if test == rare_pair:
            return hap_idx
    return None


PAIR_LABELS = {0: "00", 1: "10", 2: "01", 3: "11"}


def summarize(pattern_a_hits, sites):
    """Raw priority score per (site_index, haplotype_index)."""
    columns = [site.alleles for site in sites]
    summary = defaultdict(lambda: defaultdict(int))
    for hit in pattern_a_hits:
        for pair_name, i, j in (("ab", hit["a"], hit["b"]), ("bc", hit["b"], hit["c"])):
            counts = hit[f"{pair_name}_counts"]
            for cell_idx, count in enumerate(counts):
                if count <= CUTOFF_FREQ:
                    rare_pair = PAIR_LABELS[cell_idx]
                    hap_idx = flag_haplotype(rare_pair, i, j, columns)
                    if hap_idx is not None:
                        summary[i][hap_idx] += 1
                        summary[j][hap_idx] += 1
    return summary


def middle_correction(n_informative_sites):
    """Sites in the middle of the window have more chances to appear in a
    flagged triplet than sites at either end; index 0..n-1 -> correction."""
    corrections = {}
    for x in range(n_informative_sites):  # x is 0-based rank among informative sites
        corrections[x] = x * (n_informative_sites - 1 - x)
    return corrections


def analyze_window(sites, labels, max_span, correction_factor_override=None):
    n = len(sites)
    if n < 3:
        return None
    pattern_a_hits, n_triplets_considered = find_pattern_a(sites, max_span)
    a_count = len(pattern_a_hits)
    total_possible_triplets = math.comb(n, 3)
    proportion = correction_factor_override
    if proportion is None:
        denom = (
            n_triplets_considered if n_triplets_considered else total_possible_triplets
        )
        proportion = a_count / denom if denom else 0.0

    summary = summarize(pattern_a_hits, sites)
    max_triplets_per_site = math.comb(
        n - 1, 2
    )  # pairs of other sites a given site can be in
    corrections = middle_correction(n)

    rows = []
    for site_idx, hap_scores in summary.items():
        site = sites[site_idx]
        for hap_idx, raw_score in hap_scores.items():
            denom = (max_triplets_per_site + corrections[site_idx]) * proportion
            corrected = round(raw_score / denom, 1) if denom else 0.0
            rows.append(
                {
                    "chrom": site.chrom,
                    "pos": site.pos,
                    "sample": (
                        labels[hap_idx] if hap_idx < len(labels) else f"hap{hap_idx}"
                    ),
                    "raw_score": raw_score,
                    "corrected_score": corrected,
                }
            )
    return {
        "n_informative_sites": n,
        "a_count": a_count,
        "proportion_pattern_a": proportion,
        "rows": rows,
    }


def make_windows(sites, window_size, step):
    for start in range(0, len(sites), step):
        window = sites[start : start + window_size]
        if len(window) < 3:
            break
        yield window
        if start + window_size >= len(sites):
            break


def warn_noncontiguous_gaps(window_idx, sites, max_gap):
    """Print a warning for any pair of consecutive SNPs in the window whose
    gap exceeds max_gap bp -- large gaps make a triplet look like a double
    recombinant/error even when nothing is actually wrong."""
    if not max_gap:
        return
    for prev_site, site in zip(sites, sites[1:]):
        gap = site.pos - prev_site.pos
        if gap > max_gap:
            print(
                f"# WARNING: window {window_idx} has non-contiguous SNPs on {site.chrom}: "
                f"{prev_site.pos} to {site.pos} is a {gap} bp gap (> --max-gap {max_gap}).",
                file=sys.stderr,
            )


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("vcf", help="Input VCF file, optionally gzipped (- for stdin)")
    parser.add_argument(
        "--window",
        type=int,
        default=100,
        help="Number of informative SNP sites per window (default: 100)",
    )
    parser.add_argument(
        "--step",
        type=int,
        default=None,
        help="Step between windows in sites (default: same as --window, i.e. non-overlapping)",
    )
    parser.add_argument(
        "--max-missing",
        type=float,
        default=0.2,
        help="Exclude a site if more than this fraction of haplotypes are missing (default: 0.2)",
    )
    parser.add_argument(
        "--max-span",
        type=int,
        default=None,
        help="Exclude a triplet if the outer two sites are more than this many bp apart (guards against spurious flags from sparse, distant SNPs)",
    )
    parser.add_argument(
        "--max-gap",
        type=int,
        default=10000,
        help="Warn if two consecutive SNPs within a window are more than this many bp apart, i.e. the window has become non-contiguous (default: 10000; use 0 to disable)",
    )
    parser.add_argument(
        "--exclude-near-indel",
        type=int,
        default=0,
        help="Exclude SNPs within this many bp of an indel record (default: 0, disabled)",
    )
    parser.add_argument(
        "--allow-unphased",
        action="store_true",
        help="Use unphased heterozygous genotypes as-is (arbitrary allele order) instead of treating them as missing. Not needed for fully-homozygous VCFs, which are always collapsed to haploid regardless of phasing.",
    )
    parser.add_argument(
        "-o", "--output", default="-", help="Output file (default: stdout)"
    )
    args = parser.parse_args()

    if args.step is None:
        args.step = args.window

    sites_gen = read_vcf_sites(
        args.vcf, args.max_missing, args.exclude_near_indel, args.allow_unphased
    )
    sites = list(sites_gen)
    samples = []
    with open_maybe_gzip(args.vcf) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                break

    if not sites:
        print("# No informative biallelic SNP sites found.", file=sys.stderr)
        sys.exit(1)

    out = sys.stdout if args.output == "-" else open(args.output, "w")
    out.write(
        "#window\tchrom\tstart_pos\tend_pos\tn_informative_sites\tn_pattern_a\tproportion_pattern_a\tphasing\n"
    )

    window_idx = 0
    for window_sites in make_windows(sites, args.window, args.step):
        window_idx += 1
        warn_noncontiguous_gaps(window_idx, window_sites, args.max_gap)

        collapsed_sites, labels, collapsed = collapse_to_haploid_if_possible(
            window_sites, samples
        )
        phasing_note = "collapsed-homozygous" if collapsed else "diploid-haplotypes"

        result = analyze_window(collapsed_sites, labels, args.max_span)
        if result is None:
            continue

        start_pos = collapsed_sites[0].pos
        end_pos = collapsed_sites[-1].pos
        chrom = collapsed_sites[0].chrom
        out.write(
            f"{window_idx}\t{chrom}\t{start_pos}\t{end_pos}\t{result['n_informative_sites']}\t"
            f"{result['a_count']}\t{result['proportion_pattern_a']:.4f}\t{phasing_note}\n"
        )
        if result["rows"]:
            out.write("#\tpos\tsample\traw_score\tcorrected_score\n")
            for row in sorted(result["rows"], key=lambda r: -r["raw_score"]):
                out.write(
                    f"\t{row['pos']}\t{row['sample']}\t{row['raw_score']}\t{row['corrected_score']}\n"
                )

    if out is not sys.stdout:
        out.close()


if __name__ == "__main__":
    main()
