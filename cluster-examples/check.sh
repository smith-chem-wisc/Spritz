#!/usr/bin/env bash
# What to assert once a case has run. Each check corresponds to something that is currently a
# structural claim only - the workflow has never been executed against a bacterial reference.
#
#   ./check.sh <case> <workdir>
set -uo pipefail
CASE_NAME="${1:?usage: check.sh <case> <workdir>}"; ROOT="${2:?}"
source "$(dirname "$0")/cases/${CASE_NAME}.env"
R="$ROOT/results"; RES="$ROOT/resources"
fail=0
ok(){ printf "  PASS  %s\n" "$1"; }
no(){ printf "  FAIL  %s\n" "$1"; fail=1; }

echo "--- ${CASE_NAME} ---"

# 1. The reference actually came from Ensembl Genomes, filenames read from CHECKSUMS.
REF="${SPECIES}.${GENOME}"
[ -s "$RES/ensembl/${REF}.dna.primary_assembly.fa" ] \
  && ok "genome FASTA downloaded" || no "genome FASTA missing"
[ -s "$RES/ensembl/${REF}.pep.all.fa" ] \
  && ok "reference proteome downloaded" || no "reference proteome missing"

# 2. The codon table. This is the headline output-changing claim and the reason the hard case
#    exists: Mollicutes need table 4, everything else table 11, and getting it wrong is silent
#    except in the protein sequences.
CFG="$RES/SnpEff/snpEff.config"
if grep -q "^\s*${REF}.codonTable : ${EXPECT_CODON_TABLE}\s*$" "$CFG" 2>/dev/null; then
  ok "codon table is ${EXPECT_CODON_TABLE}"
else
  no "codon table not ${EXPECT_CODON_TABLE}; got: $(grep -F "${REF}.codonTable" "$CFG" 2>/dev/null || echo none)"
fi

# 3. The decisive one, and it needs no proteomics. SnpEff translates the gene model itself; if the
#    codon table is wrong those proteins disagree with the pep.all.fa Ensembl ships for the same
#    genome. For a Mollicute translated with table 11 every protein truncates at its first TGA, so
#    the median length collapses.
PROT=$(ls "$R"/variants/*.protein.fasta 2>/dev/null | head -1)
if [ -n "$PROT" ] && [ -s "$RES/ensembl/${REF}.pep.all.fa" ]; then
  python3 - "$PROT" "$RES/ensembl/${REF}.pep.all.fa" <<'PY'
import sys, statistics
def lens(p):
    out, cur = [], 0
    for line in open(p):
        if line.startswith(">"):
            if cur: out.append(cur)
            cur = 0
        else: cur += len(line.strip())
    if cur: out.append(cur)
    return out
a, b = lens(sys.argv[1]), lens(sys.argv[2])
ma, mb = statistics.median(a), statistics.median(b)
ratio = ma / mb if mb else 0
print(f"  ..... spritz median protein length {ma:.0f} over {len(a)} seqs; "
      f"Ensembl {mb:.0f} over {len(b)}; ratio {ratio:.2f}")
print(("  PASS  protein lengths agree with Ensembl (ratio within 0.9-1.1)" if 0.9 <= ratio <= 1.1
       else "  FAIL  protein lengths disagree - wrong codon table truncates at the first TGA"))
PY
else
  echo "  ..... no protein FASTA to compare (expected only if the run did not reach annotation)"
fi

# 4. The proteome lookup. Bacteria resolve by NCBI taxonomy id; the name-scanning path cannot
#    reach them, because its query returns only UniProt's default first page of 25 organisms.
#
#    Checking that the file EXISTS is not enough, and this check used to do only that. The rule
#    names the output after the run's species, so the path is right whatever the script downloaded:
#    the first cluster run passed this check holding 25 Homo sapiens proteins, because the scripts
#    were reading the packaged default config. Assert on the contents.
[ -s "$RES/uniprot/${SPECIES}.protein.xml.gz" ] \
  && ok "UniProt proteome downloaded" || no "UniProt proteome missing - check get_proteome"

UFA="$RES/uniprot/${SPECIES}.protein.fasta"
if [ -s "$UFA" ]; then
  uorg=$(grep -o 'OS=[^=]*' "$UFA" | sed 's/OS=//; s/ [A-Z][A-Z]*$//' \
           | sort | uniq -c | sort -rn | head -1 | sed 's/^ *[0-9]* *//')
  ucount=$(grep -c '^>' "$UFA")
  echo "  ..... UniProt proteome: ${ucount} seqs, dominant organism '${uorg}'"
  case "$uorg" in
    *"Homo sapiens"*|*"Mus musculus"*|*"Saccharomyces"*)
      no "UniProt proteome is ${uorg}, not this run's organism - scripts read the wrong config" ;;
    "") no "could not read an OS= field from ${UFA}" ;;
    *)  ok "UniProt proteome is for ${uorg}" ;;
  esac
  # 25 is UniProt's default page size, and the exact signature of the name-scanning path returning
  # its first page rather than a resolved proteome. No bacterial proteome is that small.
  [ "$ucount" -gt 25 ] && ok "UniProt proteome is not a default first page" \
    || no "UniProt proteome has ${ucount} seqs - looks like an unresolved first page"
else
  no "UniProt proteome FASTA missing at ${UFA}"
fi

# 5. The generated config records what the run actually decided.
CONF="$R/config/config.yaml"
grep -q 'division: "bacteria"' "$CONF" 2>/dev/null && ok "division recorded as bacteria" \
  || no "division not recorded as bacteria"
if [ -z "${VCF:-}" ] && [ -n "${SRA:-}" ]; then
  grep -q 'known_sites: "bootstrap"' "$CONF" 2>/dev/null \
    && ok "known_sites coerced to bootstrap" || no "known_sites not bootstrap"
  [ -s "$R/variants/bootstrap.knownsites.vcf" ] \
    && ok "bootstrap known sites produced" || no "bootstrap known sites missing"
fi

# 6. prose.txt must describe what ran. It reads the config through SPRITZ_CONFIG now; before that
#    fix it read the packaged default and described a Homo_sapiens quant run.
if [ -s "$R/prose.txt" ]; then
  grep -qi "homo sapiens\|human" "$R/prose.txt" && no "prose.txt mentions human on a bacterial run" \
    || ok "prose.txt does not claim human"
  if [ -n "${SRA:-}" ]; then
    grep -qi "bootstrapped set" "$R/prose.txt" && ok "prose.txt discloses bootstrapped recalibration" \
      || no "prose.txt omits the bootstrap disclosure"
  fi
fi

# 7. The database the user actually consumes.
[ -s "$R/final/${REF}.${RELEASE}.protein.withmods.xml.gz" ] \
  && ok "final reference database present" \
  || echo "  ..... no final reference database (only expected when the run requested one)"

echo "--- ${CASE_NAME}: $([ $fail -eq 0 ] && echo ALL CHECKS PASSED || echo "FAILURES ABOVE") ---"
exit $fail
