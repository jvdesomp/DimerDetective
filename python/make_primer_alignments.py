"""
make_primer_alignments.py
─────────────────────────
Computes primer3 nearest-neighbor 3'-end ΔG for primer pairs or a multiplex
primer set and writes alignment and/or results files.

════════════════════════════════════════════════════════════════════════════════
PAIR MODE  (--mode pairs, default)
════════════════════════════════════════════════════════════════════════════════
Reads an Excel/CSV with defined FWD/REV pairs and writes an alignment .txt.

Both heterodimer directions are computed:
  end_left-right_dg : ΔG for 3' end of FWD annealed to REV
  end_right-left_dg : ΔG for 3' end of REV annealed to FWD
  min_end_hetero     : min(end_left-right_dg, end_right-left_dg)  ← used for classification

Self-dimer end ΔG is also computed for each primer, and optionally displayed
use --type both  or  --type homodimers to include homodimer blocks.

Expected Excel columns (case-insensitive):
  forward_primer_name, forward_primer, reverse_primer_name, reverse_primer
  outcome  (optional; 'dimer' / 'clean' — used for label and prediction check)

Usage:
  python make_primer_alignments.py --input pairs.xlsx --output alignments.txt
  python make_primer_alignments.py --input pairs.xlsx --output alignments.txt \\
      --threshold -1.0 --type both

════════════════════════════════════════════════════════════════════════════════
MULTIPLEX MODE  (--mode multiplex)
════════════════════════════════════════════════════════════════════════════════
Reads an Excel/CSV with a flat primer list and computes every pairwise
combination:
  • heterodimers : every unique (i, j) pair with i < j, both directions
  • homodimers   : each primer against itself

For N primers this yields N*(N-1)/2 heterodimer pairs + N homodimers.

Writes two files:
  --output      results table (.tsv or .xlsx, detected from file extension)
  --alignments  alignment .txt  (default: <output_stem>_alignments.txt)

Table columns:
  type         heterodimer | homodimer
  id_A         primer A identifier
  seq_A        primer A sequence
  id_B         primer B identifier  (== id_A for homodimers)
  seq_B        primer B sequence    (== seq_A for homodimers)
  end_left-right_dg   ΔG for 3' end of A vs B  (kcal/mol)
  end_right-left_dg   ΔG for 3' end of B vs A  (kcal/mol; == end_left-right_dg for homodimers)
  end_left-left_dg    ΔG for 3' end of A vs itself  (kcal/mol)
  end_right-right_dg  ΔG for 3' end of B vs itself  (kcal/mol)
  min_end_hetero      min(end_left-right_dg, end_right-left_dg)
  min_end_homo        min(end_left-left_dg, end_right-right_dg)
  flagged             YES / NO  (min_end_hetero < threshold)

Expected Excel/CSV input (case-insensitive column names, all optional):
  primer_name  (or primer_id / name / id)   — omit to auto-generate P001, P002, …
  primer_seq   (or primer_sequence / sequence / seq)

  Header row is optional.  The script auto-detects whether a header is present
  and which column holds sequences (by name first, then by DNA content).

Usage:
  python make_primer_alignments.py --mode multiplex \\
      --input primers.xlsx --output results.tsv
  python make_primer_alignments.py --mode multiplex \\
      --input primers.xlsx --output results.xlsx \\
      --alignments alignments.txt --threshold -1.0 --sort asc

════════════════════════════════════════════════════════════════════════════════
Alignment display:
  Top strand    : 5'→3'
  Bottom strand : 3'→5' left to right
  | = Watson-Crick pair   . = bulge on opposite strand

  Alignment method (in order of preference):
    1. primer3-py's own thal() — the C function behind calc_end_stability,
       called with structure output.  The structure is the one behind the
       reported ΔG, with bulges and internal loops.  Where primer3 finds no
       structure (ΔG = 0) the report says [no significant structure].
    2. ntthal -a END1 (primer3 binary on PATH) — only if 1 is unavailable
       (for example on Windows).  A separately installed primer3 release may
       find another structure.  Without either, the script stops.

primer3 parameters:
  mv_conc=50 mM, dv_conc=1.5 mM, dntp_conc=0.6 mM, dna_conc=50 nM, temp_c=37 °C

Requirements:
  pip install primer3-py openpyxl pandas   (primer3-py 2.x)
"""

import argparse
import ctypes
import itertools
import re
import shutil
import subprocess
import sys
import pandas as pd
import primer3

# ntthal binary (part of the primer3 package) — fallback for alignments when
# primer3-py's thal() cannot be called directly.
_NTTHAL = shutil.which('ntthal')


# ── Thermodynamic helpers ──────────────────────────────────────────────────

PRIMER3_PARAMS = dict(mv_conc=50.0, dv_conc=1.5, dntp_conc=0.6,
                      dna_conc=50.0, temp_c=37.0)

def safe_p3(func, *args):
    try:
        return func(*args, **PRIMER3_PARAMS)
    except Exception:
        return None

def end_dg(seq1, seq2):
    """ΔG (kcal/mol) for 3' end of seq1 annealed to seq2."""
    r = safe_p3(primer3.calc_end_stability, seq1.upper(), seq2.upper())
    return r.dg / 1000 if r else 0.0


# ── Direct access to primer3-py's thal() for structures ───────────────────
# calc_end_stability() calls the C function thal() with alignment type END1
# but does not return the structure.  The compiled primer3-py module exports
# thal(), so it is called here with the same arguments plus structure output
# (as ntthal does).  The struct layouts follow thal.h in primer3-py 2.x; a
# self-test against calc_end_stability() guards against any mismatch.

class _ThalArgs(ctypes.Structure):
    _fields_ = [('type', ctypes.c_int), ('maxLoop', ctypes.c_int),
                ('mv', ctypes.c_double), ('dv', ctypes.c_double),
                ('dntp', ctypes.c_double), ('dna_conc', ctypes.c_double),
                ('temp', ctypes.c_double), ('dimer', ctypes.c_int)]

class _ThalResults(ctypes.Structure):
    _fields_ = [('msg', ctypes.c_char * 255), ('no_structure', ctypes.c_int),
                ('temp', ctypes.c_double), ('ds', ctypes.c_double),
                ('dh', ctypes.c_double), ('dg', ctypes.c_double),
                ('align_end_1', ctypes.c_int), ('align_end_2', ctypes.c_int),
                ('sec_struct', ctypes.c_void_p)]

_THAL_END1, _THL_GENERAL, _MAX_LOOP = 2, 1, 30


def _load_thal():
    try:
        from primer3 import thermoanalysis
        fn = ctypes.CDLL(thermoanalysis.__file__).thal
        fn.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.POINTER(_ThalArgs),
                       ctypes.c_int, ctypes.POINTER(_ThalResults), ctypes.c_int]
        fn.restype = None
    except Exception:
        return None
    for s1, s2 in [('CTCTATGTGTTGGAGTCATGAAG', 'CATTCCTTCTGAGGCTCATT'),
                   ('ACAACATCGTAAGGAGGAATG', 'ACAACATCGTAAGGAGGAATG')]:
        dg, _ = _call_thal(fn, s1, s2)
        if dg != primer3.calc_end_stability(s1, s2, **PRIMER3_PARAMS).dg:
            return None
    return fn


def _call_thal(fn, seq1, seq2):
    # Returns (ΔG in cal/mol, SEQ/STR structure block or '').
    p = PRIMER3_PARAMS
    args = _ThalArgs(_THAL_END1, _MAX_LOOP, p['mv_conc'], p['dv_conc'], p['dntp_conc'],
                     p['dna_conc'], p['temp_c'] + 273.15, 1)
    res = _ThalResults()
    buf = ctypes.create_string_buffer((len(seq1) + len(seq2)) * 4 + 64)
    res.sec_struct = ctypes.addressof(buf)
    fn(seq1.encode(), seq2.encode(), ctypes.byref(args), _THL_GENERAL, ctypes.byref(res), 1)
    if res.msg or res.no_structure:
        return res.dg, ''
    return res.dg, buf.value.decode()


_THAL = _load_thal()


# ── Sequence helpers ───────────────────────────────────────────────────────

def complement(base):
    return {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}.get(base.upper(), 'N')

def reverse_complement(seq):
    return ''.join(complement(b) for b in reversed(seq.upper()))


# ── Alignment ─────────────────────────────────────────────────────────────

# ── ntthal-based alignment (exact primer3 structure) ──────────────────────

def _parse_ntthal_structure(output):
    # Parse the 4-line SEQ/STR block from ntthal stdout/stderr.
    # Returns a 3-line alignment dict, or None if parsing fails.
    #
    # ntthal format (tab-separated on Linux, 4-char space prefix on macOS):
    #   SEQ<sep><line0>  top 5' overhang + bulge nucleotides + dash placeholders
    #   SEQ<sep><line1>  top strand paired blocks only
    #   STR<sep><line2>  bottom strand paired blocks only (3'->5' direction)
    #   STR<sep><line3>  bottom 3' overhang + gap markers (-) + 5' overhang
    #
    # Rules:
    #   A column is "paired" only when line1 or line2 has a non-space char.
    #   '-' in line0/line3 inside the paired range = bulge gap → rendered as '.'
    #   '|' / '.' symbols are emitted only in paired columns.
    seq_lines, str_lines = [], []
    for line in output.splitlines():
        s = line.strip()
        if s.startswith('SEQ'):
            content = line.split('\t', 1)[1] if '\t' in line else line[4:]
            seq_lines.append(content)
        elif s.startswith('STR'):
            content = line.split('\t', 1)[1] if '\t' in line else line[4:]
            str_lines.append(content)

    if len(seq_lines) != 2 or len(str_lines) != 2:
        return None

    lines = seq_lines + str_lines          # [line0, line1, line2, line3]
    W     = max(len(l) for l in lines)
    lines = [l.ljust(W) for l in lines]

    # Paired columns: where lines[1] or lines[2] carry actual bases.
    # Dashes in lines[0]/[3] are genuine bulge markers only when the column
    # sits strictly between the two outermost paired columns; otherwise they
    # are overhang placeholders and should render as a space.
    paired = [c for c in range(W) if lines[1][c] != ' ' or lines[2][c] != ' ']
    p_min  = min(paired) if paired else W
    p_max  = max(paired) if paired else -1

    top, bot = [], []
    for c in range(W):
        t0, t1 = lines[0][c], lines[1][c]
        b2, b3 = lines[2][c], lines[3][c]
        # Top: prefer paired-row (line1); '-' in line0 inside paired range = bulge
        if t1 != ' ':
            t = t1
        elif t0 == '-':
            t = '.' if p_min < c < p_max else ' '
        elif t0 != ' ':
            t = t0
        else:
            t = ' '
        # Bottom: same logic mirrored
        if b2 != ' ':
            b = b2
        elif b3 == '-':
            b = '.' if p_min < c < p_max else ' '
        elif b3 != ' ':
            b = b3
        else:
            b = ' '
        top.append(t)
        bot.append(b)

    WC = {('A', 'T'), ('T', 'A'), ('G', 'C'), ('C', 'G')}
    match = []
    for c, (t, b) in enumerate(zip(top, bot)):
        if lines[1][c] != ' ' or lines[2][c] != ' ':
            match.append('|' if (t.upper(), b.upper()) in WC else '.')
        else:
            match.append(' ')

    top_s   = ''.join(top).rstrip()
    match_s = ''.join(match).rstrip()
    bot_s   = ''.join(bot).rstrip()

    consec_3p = 0
    for ch in reversed(match_s):
        if ch == '|':
            consec_3p += 1
        else:
            break

    return {
        'top_line':   "5'-" + top_s   + "-3'",
        'match_line': "   " + match_s,
        'bot_line':   "3'-" + bot_s   + "-5'",
        'wc_pairs':   match_s.count('|'),
        'consec_3p':  consec_3p,
        'method':     'ntthal',
    }


def _build_alignment_ntthal(top, bot):
    # Call ntthal -a END1 and parse the structure. Returns None on failure.
    if _NTTHAL is None:
        return None
    p = PRIMER3_PARAMS
    cmd = [
        _NTTHAL, '-a', 'END1',
        '-s1', top.upper(), '-s2', bot.upper(),
        '-mv', str(p['mv_conc']),
        '-dv', str(p['dv_conc']),
        '-n',  str(p['dntp_conc']),
        '-d',  str(p['dna_conc']),
        '-t',  str(p['temp_c']),
    ]
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=10)
        return _parse_ntthal_structure(r.stdout + r.stderr)
    except Exception:
        return None


def build_alignment(top, bot):
    # Return a 3-line alignment dict for the 3' end of top annealed to bot,
    # or None when primer3 finds no structure.
    # Uses primer3-py's thal() (the computation behind the reported ΔG); if
    # that is unavailable, ntthal on PATH.
    if _THAL is not None:
        _, structure = _call_thal(_THAL, top.upper(), bot.upper())
        return _parse_ntthal_structure(structure) if structure else None
    return _build_alignment_ntthal(top, bot)


def alignment_method():
    if _THAL is not None:
        return "primer3-py thal() (same computation as ΔG)"
    if _NTTHAL is not None:
        return f"ntthal at {_NTTHAL} (separate primer3 program; may differ from primer3-py)"
    return None


# ── Per-pair computation ───────────────────────────────────────────────────

def compute_pair(seq_a, seq_b):
    """
    Compute end ΔG values and both direction alignments for a primer pair.

    Works for both heterodimers (seq_a != seq_b) and homodimers (seq_a == seq_b).
    For homodimers, end_left-right_dg == end_right-left_dg == end_left-left_dg.

    Returns a dict with scores and alignment data.
    """
    seq_a = seq_a.upper()
    seq_b = seq_b.upper()

    het_fwd_dg = end_dg(seq_a, seq_b)   # 3' of left vs right
    het_rev_dg = end_dg(seq_b, seq_a)   # 3' of right vs left
    min_end_hetero = min(het_fwd_dg, het_rev_dg)

    # Self-dimer terms (kept for pair-mode backward compatibility)
    end_ll_dg   = end_dg(seq_a, seq_a)
    end_rr_dg   = end_dg(seq_b, seq_b)
    min_end_homo = min(end_ll_dg, end_rr_dg)

    # Dominant direction: lower (more negative) ΔG drives the primary alignment
    if het_fwd_dg <= het_rev_dg:   # left-right is dominant
        top_dom, bot_dom = seq_a, seq_b
        top_alt, bot_alt = seq_b, seq_a
        which = 'A'             # 3' end of A is dominant
    else:
        top_dom, bot_dom = seq_b, seq_a
        top_alt, bot_alt = seq_a, seq_b
        which = 'B'             # 3' end of B is dominant

    aln_dom = build_alignment(top_dom, bot_dom)
    # For homodimers the two directions are identical; skip the redundant call
    aln_alt = None if seq_a == seq_b else build_alignment(top_alt, bot_alt)

    return {
        'end_left-right_dg': het_fwd_dg,
        'end_right-left_dg': het_rev_dg,
        'min_end_hetero':     min_end_hetero,
        'end_left-left_dg':   end_ll_dg,
        'end_right-right_dg': end_rr_dg,
        'min_end_homo':       min_end_homo,
        'which':       which,      # 'A' or 'B' — which primer's 3' end is dominant
        'aln':         aln_dom,    # dominant direction alignment
        'aln_alt':     aln_alt,    # non-dominant direction alignment (None for homodimers)
    }


# ── Text-output helpers ────────────────────────────────────────────────────

_SEP = '─' * 72

def _aln_block(aln, top_id, bot_id, indent='  '):
    """
    Return indented alignment lines (list) labeled with primer IDs.
    top_id : identifier of the primer whose 3' end is shown on top.
    bot_id : identifier of the primer on the bottom strand.
    """
    if not aln:
        return [f"{indent}alignment of 3' {top_id} vs {bot_id}: [no significant structure]"]
    return [
        f"{indent}alignment of 3' {top_id} vs {bot_id} "
        f"({aln['wc_pairs']} WC pairs, "
        f"{aln['consec_3p']} consecutive at 3' end)",
        f"{indent}{aln['top_line']}",
        f"{indent}{aln['match_line']}",
        f"{indent}{aln['bot_line']}",
    ]


# ── Pair mode formatter ────────────────────────────────────────────────────

def format_pair(row, result, threshold, show_homodimers=False):
    """Format one entry for pair-mode alignment .txt output."""
    fwd_name = row['fwd_name']
    rev_name = row['rev_name']
    lines    = [_SEP]

    lines.append(f"{fwd_name} / {rev_name}  [heterodimer]")
    lines.append('')

    call = 'dimer' if result['min_end_hetero'] < threshold else 'clean'
    if 'outcome' in row and pd.notna(row['outcome']):
        label   = row['outcome'].lower()
        correct = '✓' if call == label else '✗ (label mismatch)'
        lines.append(f"  predicted: {call}   label: {label}   {correct}")
    else:
        lines.append(f"  predicted: {call}")

    dominant_id  = fwd_name if result['which'] == 'A' else rev_name
    lines.append(f"  end_left-right_dg = {result['end_left-right_dg']:7.3f} kcal/mol  (3' of {fwd_name} vs {rev_name})")
    lines.append(f"  end_right-left_dg = {result['end_right-left_dg']:7.3f} kcal/mol  (3' of {rev_name} vs {fwd_name})")
    lines.append(f"  min_end_hetero = {result['min_end_hetero']:7.3f} kcal/mol  ← dominant ({dominant_id} 3'-end)")
    lines.append('')
    lines.append(f"  {fwd_name}: 5'-{row['fwd']}-3'")
    lines.append(f"  {rev_name}: 5'-{row['rev']}-3'")
    lines.append('')

    # Dominant direction
    if result['which'] == 'A':
        dom_top, dom_bot = fwd_name, rev_name
        alt_top, alt_bot = rev_name, fwd_name
    else:
        dom_top, dom_bot = rev_name, fwd_name
        alt_top, alt_bot = fwd_name, rev_name

    lines.extend(_aln_block(result['aln'],     dom_top, dom_bot))
    lines.append('')
    lines.extend(_aln_block(result['aln_alt'], alt_top, alt_bot))
    lines.append('')

    if show_homodimers:
        lines.append(f"  end_left-left_dg ({fwd_name}): {result['end_left-left_dg']:.3f} kcal/mol")
        lines.append(f"  end_right-right_dg ({rev_name}): {result['end_right-right_dg']:.3f} kcal/mol")
        lines.append('')

    return lines


def format_homodimer_entry(pid, seq, result):
    """Format a homodimer block for pair mode --homodimers output."""
    lines = [_SEP]
    lines.append(f"{pid} / {pid}  [homodimer]")
    lines.append('')
    lines.append(f"  min_end_homo = {result['min_end_hetero']:7.3f} kcal/mol")
    lines.append('')
    lines.append(f"  {pid}: 5'-{seq}-3'")
    lines.append('')
    lines.extend(_aln_block(result['aln'], pid, pid))
    lines.append('')
    return lines


# ── Multiplex mode formatter ───────────────────────────────────────────────

def format_multiplex_entry(pair_type, id_a, seq_a, id_b, seq_b,
                           result, threshold):
    """Format one entry for multiplex alignment .txt output."""
    lines = [_SEP]

    if pair_type == 'homodimer':
        lines.append(f"{id_a} / {id_a}  [homodimer]")
        lines.append('')
        lines.append(f"  predicted: {'dimer' if result['min_end_hetero'] < threshold else 'clean'}")
        lines.append(f"  min_end_homo = {result['min_end_hetero']:7.3f} kcal/mol")
        lines.append('')
        lines.append(f"  {id_a}: 5'-{seq_a}-3'")
        lines.append('')
        lines.extend(_aln_block(result['aln'], id_a, id_a))
    else:
        dominant_id  = id_a if result['which'] == 'A' else id_b
        lines.append(f"{id_a} / {id_b}  [heterodimer]")
        lines.append('')
        lines.append(f"  predicted: {'dimer' if result['min_end_hetero'] < threshold else 'clean'}")
        lines.append(f"  end_left-right_dg = {result['end_left-right_dg']:7.3f} kcal/mol  (3' of {id_a} vs {id_b})")
        lines.append(f"  end_right-left_dg = {result['end_right-left_dg']:7.3f} kcal/mol  (3' of {id_b} vs {id_a})")
        lines.append(f"  min_end_hetero    = {result['min_end_hetero']:7.3f} kcal/mol  ← dominant ({dominant_id} 3'-end)")
        lines.append(f"  end_left-left_dg  = {result['end_left-left_dg']:7.3f} kcal/mol  (3' of {id_a} vs {id_a})")
        lines.append(f"  end_right-right_dg= {result['end_right-right_dg']:7.3f} kcal/mol  (3' of {id_b} vs {id_b})")
        lines.append(f"  min_end_homo      = {result['min_end_homo']:7.3f} kcal/mol")
        lines.append('')
        lines.append(f"  {id_a}: 5'-{seq_a}-3'")
        lines.append(f"  {id_b}: 5'-{seq_b}-3'")
        lines.append('')

        # Dominant alignment first, then alternative
        if result['which'] == 'A':
            dom_top, dom_bot = id_a, id_b
            alt_top, alt_bot = id_b, id_a
        else:
            dom_top, dom_bot = id_b, id_a
            alt_top, alt_bot = id_a, id_b

        lines.extend(_aln_block(result['aln'],     dom_top, dom_bot))
        lines.append('')
        lines.extend(_aln_block(result['aln_alt'], alt_top, alt_bot))

    lines.append('')
    return lines


# ── Column name normalization ──────────────────────────────────────────────

def find_col(df, candidates):
    """Return the first matching column name (case-insensitive), or None.
    Tolerates integer column indices produced by header=None reads."""
    lower = {str(c).lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


# ── Flexible primer-list reader ────────────────────────────────────────────

_DNA_PAT = re.compile(r'^[ATGCUNatgcun\-]+$')

def _looks_like_dna(val, min_len=6):
    """Return True if val is a plausible DNA/RNA sequence."""
    s = str(val).strip()
    return len(s) >= min_len and bool(_DNA_PAT.match(s))

def _dna_fraction(series, min_len=6):
    """Fraction of non-null values in series that look like DNA sequences."""
    vals = series.dropna().astype(str).str.strip()
    return sum(_looks_like_dna(v, min_len) for v in vals) / len(vals) if len(vals) else 0.0

def _read_primer_list(path):
    """
    Read a flat primer list from an Excel or CSV file.  Handles four layouts:

      1. Header + ID column + seq column  → columns matched by name
      2. Header + seq column only         → IDs auto-generated as P001, P002, …
      3. No header + 1 column             → seq detected by DNA content; IDs auto-generated
      4. No header + 2 columns            → seq detected by DNA content; other col = ID

    Returns a list of (primer_id, primer_seq) tuples.
    """
    read_fn = pd.read_csv if path.lower().endswith('.csv') else pd.read_excel

    df_peek = read_fn(path)

    _label_candidates = ['primer_seq', 'primer_sequence', 'sequence', 'seq',
                         'primer_name', 'primer_id', 'name', 'id']
    any_recognized = find_col(df_peek, _label_candidates) is not None
    all_dna        = all(_looks_like_dna(str(c)) for c in df_peek.columns)
    any_dna        = any(_looks_like_dna(str(c)) for c in df_peek.columns)

    has_header = any_recognized and not all_dna
    if not has_header and not any_dna:
        has_header = True

    df = df_peek if has_header else read_fn(path, header=None)

    # ── Locate the sequence column ─────────────────────────────────────────
    seq_col = find_col(df, ['primer_seq', 'primer_sequence', 'sequence', 'seq'])
    if seq_col is None:
        scores = {c: _dna_fraction(df[c]) for c in df.columns}
        best   = max(scores, key=scores.get)
        if scores[best] < 0.8:
            sys.exit(
                f"ERROR: could not identify a primer sequence column.\n"
                f"DNA-content scores per column: {scores}\n"
                f"Tip: add a header row with a column named 'primer_seq'."
            )
        seq_col = best

    # ── Locate the ID column ───────────────────────────────────────────────
    id_col = find_col(df, ['primer_name', 'primer_id', 'name', 'id', 'ID'])
    if id_col is None:
        others  = [c for c in df.columns if c != seq_col]
        non_dna = [c for c in others if _dna_fraction(df[c]) < 0.5]
        id_col  = non_dna[0] if non_dna else None

    # ── Build result list ──────────────────────────────────────────────────
    n_rows   = len(df)
    n_digits = len(str(n_rows))
    primers  = []
    for _, row in df.iterrows():
        seq = str(row[seq_col]).strip().upper()
        if not seq or seq in ('NAN', 'NONE'):
            continue
        if id_col is not None:
            raw = row[id_col]
            try:
                pid = str(int(raw))
            except (ValueError, TypeError):
                pid = str(raw).strip()
            if pid.lower() in ('nan', 'none', ''):
                pid = f"P{len(primers) + 1:0{n_digits}d}"
        else:
            pid = f"P{len(primers) + 1:0{n_digits}d}"
        primers.append((pid, seq))

    return primers


# ── Mode: pairs ────────────────────────────────────────────────────────────

def run_pairs(args):
    if args.input.lower().endswith('.csv'):
        df = pd.read_csv(args.input)
    else:
        df = pd.read_excel(args.input)

    col_fwd_name = find_col(df, [
        'forward_primer_name', 'fwd_name', 'fwd_id', 'forward_name', 'forward_id',
        'left_primer_name',    'left_name', 'left_id', 'left-id',
        'Forward Primer Name',
    ])
    col_fwd_seq  = find_col(df, [
        'forward_primer', 'forward_primer_seq', 'forward_primer_sequence',
        'forward_sequence', 'fwd', 'fwd_seq',
        'left_primer',    'left_primer_seq',    'left_primer_sequence',
        'left_sequence',  'left-sequence', 'left-primer',
        'Forward Primer Seq', 'Forward Primer',
    ])
    col_rev_name = find_col(df, [
        'reverse_primer_name', 'rev_name', 'rev_id', 'reverse_name', 'reverse_id',
        'right_primer_name',   'right_name', 'right_id', 'right-id',
        'Reverse Primer Name',
    ])
    col_rev_seq  = find_col(df, [
        'reverse_primer', 'reverse_primer_seq', 'reverse_primer_sequence',
        'reverse_sequence', 'rev', 'rev_seq',
        'right_primer',   'right_primer_seq',   'right_primer_sequence',
        'right_sequence', 'right-sequence', 'right-primer',
        'Reverse Primer Seq', 'Reverse Primer',
    ])
    col_outcome  = find_col(df, ['outcome', 'label', 'class', 'result', 'status'])

    # ── DNA-content fallback when column names are not recognized ──────────
    # If either sequence column was not found by name, try to identify the
    # two sequence columns by DNA content, then assign left→fwd, right→rev.
    if col_fwd_seq is None or col_rev_seq is None:
        scores   = {c: _dna_fraction(df[c]) for c in df.columns}
        dna_cols = sorted([c for c, s in scores.items() if s >= 0.8],
                          key=lambda c: list(df.columns).index(c))
        if len(dna_cols) >= 2:
            col_fwd_seq = col_fwd_seq or dna_cols[0]
            col_rev_seq = col_rev_seq or dna_cols[1]
            # ID columns: non-DNA columns adjacent to the sequence columns
            non_dna = [c for c in df.columns if scores[c] < 0.5
                       and c not in (col_fwd_seq, col_rev_seq)
                       and (col_outcome is None or c != col_outcome)]
            if col_fwd_name is None and len(non_dna) >= 1:
                col_fwd_name = non_dna[0]
            if col_rev_name is None and len(non_dna) >= 2:
                col_rev_name = non_dna[1]

    for col, label in [(col_fwd_seq, 'forward/left primer sequence'),
                       (col_rev_seq,  'reverse/right primer sequence')]:
        if col is None:
            sys.exit(
                f"ERROR: could not find column for {label}.\n"
                f"Available columns: {list(df.columns)}\n"
                f"Recognized name patterns: forward/fwd/left + _primer/_sequence/_seq\n"
                f"                          reverse/rev/right + _primer/_sequence/_seq\n"
                f"If your columns have other names, rename them or add a row with these headers."
            )

    # ── Derive alignment output path if not given ────────────────────────
    align_path = args.alignments
    if align_path is None:
        stem       = args.output.rsplit('.', 1)[0] if '.' in args.output else args.output
        align_path = stem + '_alignments.txt'

    # ── Compute all pairs ──────────────────────────────────────────────────
    n_pairs       = len(df)
    n_flagged_het = 0
    pair_rows     = []
    tsv_rows      = []

    for i, df_row in df.iterrows():
        fwd      = str(df_row[col_fwd_seq]).strip()
        rev      = str(df_row[col_rev_seq]).strip()
        fwd_name = str(df_row[col_fwd_name]).strip() if col_fwd_name else f"FWD_{i + 1}"
        rev_name = str(df_row[col_rev_name]).strip() if col_rev_name else f"REV_{i + 1}"
        outcome  = df_row[col_outcome] if col_outcome else None
        row = {
            'fwd_name': fwd_name,
            'rev_name': rev_name,
            'fwd':      fwd,
            'rev':      rev,
            'outcome':  outcome,
        }
        result  = compute_pair(fwd, rev)
        flagged = result['min_end_hetero'] < args.threshold
        if flagged:
            n_flagged_het += 1
        pair_rows.append((row, result))

        tsv_row = {
            'fwd_name':           fwd_name,
            'fwd_seq':            fwd,
            'rev_name':           rev_name,
            'rev_seq':            rev,
            'end_left-right_dg':  round(result['end_left-right_dg'],  4),
            'end_right-left_dg':  round(result['end_right-left_dg'],  4),
            'min_end_hetero':     round(result['min_end_hetero'],      4),
            'end_left-left_dg':   round(result['end_left-left_dg'],   4),
            'end_right-right_dg': round(result['end_right-right_dg'], 4),
            'min_end_homo':       round(result['min_end_homo'],        4),
            'flagged':            'YES' if flagged else 'NO',
            'predicted':          'dimer' if flagged else 'clean',
        }
        if col_outcome:
            tsv_row['outcome'] = outcome
        tsv_rows.append(tsv_row)

    # ── Write results table ────────────────────────────────────────────────
    out_df = pd.DataFrame(tsv_rows)
    if args.output.lower().endswith('.xlsx'):
        out_df.to_excel(args.output, index=False)
        print(f"Written {n_pairs} pairs → {args.output}  (Excel)")
    else:
        out_df.to_csv(args.output, sep='\t', index=False)
        print(f"Written {n_pairs} pairs → {args.output}  (TSV)")

    # ── Write alignment text ───────────────────────────────────────────────
    pair_type = getattr(args, 'type', 'both')
    aln_lines = [
        '=' * 72,
        'DimerDetective alignment output',
        '',
        f"mode      : pairs",
        f"type      : {pair_type}",
        f"threshold : {args.threshold:.2f} kcal/mol → dimer",
        f"pairs     : {n_pairs}   |   flagged: {n_flagged_het}",
        "| = Watson-Crick pair",
        ". = bulge on opposite strand",
        '=' * 72,
        '',
    ]

    if pair_type in ('heterodimers', 'both'):
        for row, result in pair_rows:
            aln_lines.extend(format_pair(row, result, args.threshold,
                                         show_homodimers=False))

    if pair_type in ('homodimers', 'both'):
        aln_lines += ['', '=' * 72, 'homodimers', '=' * 72, '']
        seen = set()
        for row, result in pair_rows:
            for pid, seq in [
                (row['fwd_name'], row['fwd']),
                (row['rev_name'], row['rev']),
            ]:
                if pid not in seen:
                    seen.add(pid)
                    hom = compute_pair(seq, seq)
                    aln_lines.extend(format_homodimer_entry(pid, seq, hom))

    with open(align_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(aln_lines))
    print(f"Written alignments → {align_path}")


# ── Mode: multiplex ────────────────────────────────────────────────────────

def _sort_label(sort):
    return {
        'asc':   'ascending (most stable dimer first)',
        'desc':  'descending (least stable dimer first)',
        'input': 'input order',
    }.get(sort, sort)


def run_multiplex(args):
    primers = _read_primer_list(args.input)

    n          = len(primers)
    n_hetero   = n * (n - 1) // 2
    mplex_type = getattr(args, 'type', 'heterodimers')

    if mplex_type == 'heterodimers':
        summary = f"{n_hetero} heterodimers = {n_hetero} computations"
    elif mplex_type == 'homodimers':
        summary = f"{n} homodimers = {n} computations"
    else:  # both
        n_total = n_hetero + n
        summary = f"{n_hetero} heterodimers + {n} homodimers = {n_total} computations"

    print(f"Primers: {n}  →  {summary}", flush=True)

    # Derive alignment output path if not given
    align_path = args.alignments
    if align_path is None:
        stem       = args.output.rsplit('.', 1)[0] if '.' in args.output else args.output
        align_path = stem + '_alignments.txt'

    # ── Compute all entries ────────────────────────────────────────────────
    entries = []   # (pair_type, id_a, seq_a, id_b, seq_b, result)

    mplex_type = getattr(args, 'type', 'both')
    if mplex_type in ('heterodimers', 'both'):
        for (id_a, seq_a), (id_b, seq_b) in itertools.combinations(primers, 2):
            result = compute_pair(seq_a, seq_b)
            entries.append(('heterodimer', id_a, seq_a, id_b, seq_b, result))

    if mplex_type in ('homodimers', 'both'):
        for pid, seq in primers:
            result = compute_pair(seq, seq)
            entries.append(('homodimer', pid, seq, pid, seq, result))

    # ── Sort ──────────────────────────────────────────────────────────────
    sort = getattr(args, 'sort', 'asc')
    if sort == 'asc':
        entries.sort(key=lambda e: e[5]['min_end_hetero'])
    elif sort == 'desc':
        entries.sort(key=lambda e: e[5]['min_end_hetero'], reverse=True)
    # 'input': preserve natural order (heterodimers first, then homodimers)

    # ── Write results table (.tsv or .xlsx) ───────────────────────────────
    tsv_rows = []
    for pair_type, id_a, seq_a, id_b, seq_b, result in entries:
        flagged = 'YES' if result['min_end_hetero'] < args.threshold else 'NO'
        tsv_rows.append({
            'type':               pair_type,
            'id_A':               id_a,
            'seq_A':              seq_a,
            'id_B':               id_b,
            'seq_B':              seq_b,
            'end_left-right_dg':  round(result['end_left-right_dg'],  4),
            'end_right-left_dg':  round(result['end_right-left_dg'],  4),
            'min_end_hetero':     round(result['min_end_hetero'],      4),
            'end_left-left_dg':   round(result['end_left-left_dg'],   4),
            'end_right-right_dg': round(result['end_right-right_dg'], 4),
            'min_end_homo':       round(result['min_end_homo'],        4),
            'flagged':            flagged,
        })

    n_flagged = sum(1 for r in tsv_rows if r['flagged'] == 'YES')
    out_df    = pd.DataFrame(tsv_rows)

    if args.output.lower().endswith('.xlsx'):
        out_df.to_excel(args.output, index=False)
        print(f"Written {len(tsv_rows)} rows → {args.output}  (Excel)")
    else:
        out_df.to_csv(args.output, sep='\t', index=False)
        print(f"Written {len(tsv_rows)} rows → {args.output}  (TSV)")

    # ── Write alignment .txt ───────────────────────────────────────────────
    aln_lines = [
        '=' * 72,
        'DimerDetective alignment output',
        '',
        f"mode      : multiplex",
        f"threshold : {args.threshold:.2f} kcal/mol → dimer",
        f"primers   : {n}   |   pairs: {len(entries)}   |   flagged: {n_flagged}",
        f"type      : {getattr(args, 'type', 'both')}",
        f"sorted by : {_sort_label(sort)}",
        "| = Watson-Crick pair",
        ". = bulge on opposite strand",
        '=' * 72,
        '',
    ]

    for pair_type, id_a, seq_a, id_b, seq_b, result in entries:
        aln_lines.extend(
            format_multiplex_entry(pair_type, id_a, seq_a, id_b, seq_b,
                                   result, args.threshold)
        )

    with open(align_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(aln_lines))
    print(f"Written alignments → {align_path}")


# ── Main ───────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        '--input', required=True,
        help='Input Excel (.xlsx) or CSV file',
    )
    parser.add_argument(
        '--output', required=True,
        help=(
            'pairs mode    : alignment .txt output file\n'
            'multiplex mode: results table; .tsv → tab-separated, .xlsx → Excel'
        ),
    )
    parser.add_argument(
        '--alignments', default=None,
        help=(
            'multiplex mode only: alignment .txt output file '
            '(default: <output_stem>_alignments.txt)'
        ),
    )
    parser.add_argument(
        '--threshold', type=float, default=-1.0,
        help='min_end_hetero threshold for DIMER flag (default: -1.0 kcal/mol)',
    )
    parser.add_argument(
        '--mode', choices=['pairs', 'multiplex'], default='pairs',
        help='Analysis mode (default: pairs)',
    )
    parser.add_argument(
        '--sort', choices=['asc', 'desc', 'input'], default='asc',
        help=(
            'multiplex mode: sort order for output\n'
            '  asc   — ascending min_end_hetero (most stable dimer first, default)\n'
            '  desc  — descending min_end_hetero (least stable first)\n'
            '  input — preserve input primer order'
        ),
    )
    parser.add_argument(
        '--type', choices=['heterodimers', 'homodimers', 'both'], default='heterodimers',
        help=(
            'which interactions to compute and report (default: both)\n'
            '  heterodimers — only cross-primer interactions\n'
            '  homodimers   — only self-dimer interactions\n'
            '  both         — heterodimers and homodimers'
        ),
    )
    args = parser.parse_args()
    method = alignment_method()
    if method is None:
        sys.exit("ERROR: cannot draw primer3 alignments: primer3-py's thal() is not\n"
                 "accessible on this system and ntthal is not on PATH.\n"
                 "Install primer3 (conda install -c bioconda primer3) to provide ntthal.")
    print(f"Alignments: {method}")

    if args.mode == 'pairs':
        run_pairs(args)
    else:
        run_multiplex(args)


if __name__ == '__main__':
    main()
