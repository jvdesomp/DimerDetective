#!/usr/bin/env python3
"""
Check that the web version gives byte-identical output to the Python script.

For each input file, runs python/make_primer_alignments.py and
tests/dd-cli.js (the web core in Node) with the same settings and compares
the results tables and the alignment reports. The check stops if the
script cannot draw its alignments with primer3-py's own thal().

Usage:
  python tests/validate.py --mode pairs     data/pairs1.xlsx data/pairs2.xlsx
  python tests/validate.py --mode multiplex data/primers.xlsx

Requirements: pip install primer3-py==2.4.0 openpyxl pandas; node on PATH;
web/thal.js built with wasm/build.sh.
"""
import argparse
import difflib
import os
import subprocess
import sys
import tempfile

import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPT = os.path.join(ROOT, 'python', 'make_primer_alignments.py')
CLI = os.path.join(ROOT, 'tests', 'dd-cli.js')


def run(cmd):
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        sys.exit(f"command failed: {' '.join(cmd)}\n{r.stderr}")
    return r.stdout


def compare(label, a_path, b_path):
    with open(a_path, encoding='utf-8') as f:
        a = f.read().splitlines()
    with open(b_path, encoding='utf-8') as f:
        b = f.read().splitlines()
    if a == b:
        print(f"  {label}: identical ({len(a)} lines)")
        return True
    # diff only a window from the first difference: a full difflib run on
    # large reports (over a million lines) takes hours
    first = next((i for i, (x, y) in enumerate(zip(a, b)) if x != y), min(len(a), len(b)))
    n_diff = sum(x != y for x, y in zip(a, b)) + abs(len(a) - len(b))
    lo = max(0, first - 2)
    diff = list(difflib.unified_diff(a[lo:first + 40], b[lo:first + 40], 'python', 'web', lineterm='', n=1))
    print(f"  {label}: differs ({len(a)} lines in python, {len(b)} in web; first at line {first + 1}, "
          f"{n_diff} positions differ)")
    print('\n'.join('    ' + d for d in diff[:40]))
    return False


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--mode', choices=['pairs', 'multiplex'], default='pairs')
    ap.add_argument('--type', choices=['heterodimers', 'homodimers', 'both'], default='both')
    ap.add_argument('--sort', choices=['asc', 'desc', 'input'], default='asc')
    ap.add_argument('--threshold', default='-1.00')
    ap.add_argument('inputs', nargs='+')
    args = ap.parse_args()

    ok = True
    with tempfile.TemporaryDirectory() as tmp:
        for path in args.inputs:
            name = os.path.splitext(os.path.basename(path))[0]
            print(f"{name} ({args.mode}, {args.type}, sort {args.sort})")
            read = pd.read_csv if path.lower().endswith('.csv') else pd.read_excel
            tsv_in = os.path.join(tmp, name + '_input.tsv')
            read(path).to_csv(tsv_in, sep='\t', index=False)

            common = ['--mode', args.mode, '--threshold', args.threshold,
                      '--sort', args.sort, '--type', args.type]
            py_out = os.path.join(tmp, name + '_py.tsv')
            js_out = os.path.join(tmp, name + '_js.tsv')
            log = run([sys.executable, SCRIPT, '--input', path, '--output', py_out] + common)
            if 'Alignments: primer3-py thal()' not in log:
                sys.exit('the script could not call primer3-py thal(): ' + log.splitlines()[0])
            run(['node', CLI, '--input', tsv_in, '--output', js_out] + common)

            ok &= compare('results table', py_out, js_out)
            ok &= compare('alignments', py_out[:-4] + '_alignments.txt', js_out[:-4] + '_alignments.txt')
    print('ALL IDENTICAL' if ok else 'DIFFERENCES FOUND')
    sys.exit(0 if ok else 1)


if __name__ == '__main__':
    main()
