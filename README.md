# DimerDetective

DimerDetective predicts primer dimer formation from the primer3 3′ end stability (`calc_end_stability`) of every primer combination. It comes as a web app that runs entirely in the browser and as a Python command-line script. Both give identical results.

Web app: https://jvdesomp.github.io/DimerDetective/

## How it works

The thermodynamics are not reimplemented. The C source of primer3's `thal` algorithm (`thal.c`) and its parameter files are taken unchanged from [primer3-py](https://github.com/libnano/primer3-py) 2.4.0 and compiled to WebAssembly with Emscripten:

- `wasm/dd_thal.c` makes the same `thal()` call as `primer3.calc_end_stability()`: alignment type `END1`, mode `THL_GENERAL`, 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer, 37 °C and a maximum loop of 30 nt. It also asks `thal()` for the secondary structure, so each alignment diagram is the structure behind the reported ΔG.
- `web/dd-core.js` ports the input parsing, scoring metrics and output formats of `python/make_primer_alignments.py`.
- `web/app.js` is the interface. Sequences never leave the browser.

## Validation

The web app was compared with `primer3.calc_end_stability` (primer3-py 2.4.0) on more than 36,000 primer combinations, including the 200 pairs of the Biogazelle training data and random sequences. The largest ΔG difference was 5 × 10⁻⁷ cal/mol, which is floating-point rounding.

`tests/validate.py` runs the web core and the Python script on the same input and compares their results tables and alignment reports. On the Biogazelle training data, in pair and multiplex mode, both are byte-identical. GitHub Actions repeats the comparison on every push, using the public example files in `tests/`.

The two implementations differ in one input case, on purpose. If multiplex mode receives a table of primer pairs (two sequence columns), the web app uses all unique primers. The script uses only the first sequence column.

To run the checks yourself:

```bash
bash wasm/build.sh             # needs Emscripten (emcc)
node tests/smoke.js            # primer3-py reference values
python tests/validate.py --mode pairs tests/example_pairs.csv
python tests/validate.py --mode multiplex tests/example_primers.csv
```

## Repository layout

```
web/            the web app (index.html, app.js, dd-core.js, style.css, img/)
wasm/           dd_thal.c wrapper, build.sh and primer3/ (thal.c, thal.h, primer3_config/ from primer3-py 2.4.0)
python/         command-line script (make_primer_alignments.py) and wrapper (dimer-prediction.sh)
tests/          smoke test, Node command line (dd-cli.js), validation against the script, example inputs
.github/        workflow that builds the WebAssembly module and publishes web/ on GitHub Pages
```

## Running locally

```bash
bash wasm/build.sh
python3 -m http.server 8000 --directory web
```

Then open http://localhost:8000.

## Deployment

Every push to `main` triggers the GitHub Actions workflow. It compiles `thal.c` with Emscripten 6.0.13, checks the result against primer3-py values and against the Python script, and publishes `web/` to GitHub Pages. One-time setup: in the repository settings, go to **Pages** and set **Source** to **GitHub Actions**.

## Python script

ΔG comes from `primer3.calc_end_stability`. For the alignment diagrams, the script calls the C function behind it, `thal()`, inside the installed primer3-py library, asking it for the structure as well. Each diagram is therefore the structure behind the reported ΔG, and the separate `ntthal` program is no longer needed. The script prints which alignment method it uses. If `thal()` cannot be reached (for example, on Windows), it falls back to `ntthal` on PATH, which may be a different primer3 release. With neither available, it stops with an error. Use primer3-py 2.4.0 for results identical to the web app.

```bash
pip install primer3-py==2.4.0 openpyxl pandas
python python/make_primer_alignments.py --input pairs.xlsx --output results.tsv --threshold -1.00 --type both
python python/make_primer_alignments.py --mode multiplex --input primers.xlsx --output results.xlsx --sort asc
```

## License and citation

GPL-2.0, as primer3 and primer3-py. © Ghent University 2026.

- Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG (2012) Primer3, new capabilities and interfaces. *Nucleic Acids Res* 40:e115.
- SantaLucia J (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95:1460–1465.
