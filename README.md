# DimerDetective

DimerDetective predicts primer dimer formation from the primer3 3′ end stability (`calc_end_stability`) of every primer combination. It comes as a web app that runs entirely in the browser and as a Python command-line script. Both give identical results.

Web app: https://jvdesomp.github.io/DimerDetective/

## How it works

The thermodynamics are not reimplemented. The C source of primer3's `thal` algorithm (`thal.c`) and its parameter files are taken unchanged from [primer3-py](https://github.com/libnano/primer3-py) 2.4.0 and compiled to WebAssembly with Emscripten:

- `wasm/dd_thal.c` makes the same `thal()` call as `primer3.calc_end_stability()`: alignment type `END1`, mode `THL_GENERAL`, 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer, 37 °C and a maximum loop of 30 nt. It also asks `thal()` for the secondary structure, so each alignment diagram is the structure behind the reported ΔG.
- `web/dd-core.js` ports the input parsing, scoring metrics and output formats of `python/make_primer_alignments.py`.
- `web/app.js` is the interface. Sequences never leave the browser.

## Usage statistics

The web app counts anonymous usage with [GoatCounter](https://www.goatcounter.com) (site code `dimerdetective`): page views, runs per mode and interaction type, the number of pairs or primers per run, and downloads. GoatCounter uses no cookies and stores no personal data. Sequences, primer names and file names are never sent. Nothing is counted when the page is opened from `localhost` or as a local file.

## Contact form

The info tab has a contact form (email address, subject, message; no file upload). Messages are relayed by [Web3Forms](https://web3forms.com) to the developer by email, with the sender's address as reply-to. They contain what the sender enters, plus the sender's IP address, which Web3Forms adds. A hidden field stops simple spam bots. A (quick) reply cannot be guaranteed.

## Homodimers: not validated

Classification uses `min_end_hetero` only. Homodimer values (`end_left-left_dg`, `end_right-right_dg`, `min_end_homo` and homodimer rows) are reported for completeness, but they are not validated. In the Biogazelle training data, the no-template control (NTC) signal cannot tell homodimers from heterodimers, and predicted homodimers had no predictive value:

- `min_end_homo` gave an AUC of 0.60, against 0.94 for `min_end_hetero`, and added nothing to `min_end_hetero` in a logistic model (p = 0.30).
- 7 of the 100 clean pairs had a strong predicted homodimer (below −2.44 kcal/mol, down to −6.46 kcal/mol) without any NTC signal.

Interpret homodimer values with caution.

## Validation

The web app was compared with `primer3.calc_end_stability` (primer3-py 2.4.0) on more than 36,000 primer combinations, including the 200 pairs of the Biogazelle training data and random sequences. The largest ΔG difference was 5 × 10⁻⁷ cal/mol, which is floating-point rounding.

`tests/validate.py` runs the web core and the Python script on the same input and compares their results tables and alignment reports. On the Biogazelle training data, in pair and multiplex mode, both are byte-identical. GitHub Actions repeats the comparison on every push, using the public example files in `tests/`.

In multiplex mode, both accept a table of primer pairs (all primers of both sequence columns are used) and analyze identical sequences once, keeping the first ID. Both report this in a note, on screen and in the alignment report.

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

Every push to `main` triggers the GitHub Actions workflow. It compiles `thal.c` with Emscripten 6.0.13, checks the result against primer3-py values and against the Python script, stamps the links to the script and style files with the commit (so browsers load each new version instead of a cached one), and publishes `web/` to GitHub Pages. One-time setup: in the repository settings, go to **Pages** and set **Source** to **GitHub Actions**.

## Python script

ΔG comes from `primer3.calc_end_stability`. For the alignment diagrams, the script calls the C function behind it, `thal()`, inside the installed primer3-py library, asking it for the structure as well. Each diagram is therefore the structure behind the reported ΔG, and the separate `ntthal` program is no longer needed. The script prints which alignment method it uses. If `thal()` cannot be reached (for example, on Windows), it falls back to `ntthal` on PATH, which may be a different primer3 release. With neither available, it stops with an error. Use primer3-py 2.4.0 for results identical to the web app.

```bash
pip install primer3-py==2.4.0 openpyxl pandas
python python/make_primer_alignments.py --input pairs.xlsx --output results.tsv --threshold -1.00 --type both
python python/make_primer_alignments.py --mode multiplex --input primers.xlsx --output results.xlsx --sort asc
```

## How to cite

DimerDetective version 1.2, Jo Vandesompele, Ghent University, https://jvdesomp.github.io/DimerDetective/

Every output file of the web app and the script states the version and build, for example `1.2 (build 20261102-0058)`: in the header of the alignment report and in a final `version` column of the tables. The build is the date and time (Belgian time, YYYYMMDD-hhmm) of the last change; small changes keep version 1.2 and get a new build. The web app shows the build in its footer; the script prints it, also with `--version`.

## License and references

GPL-2.0, as primer3 and primer3-py. © Ghent University 2026.

- Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG (2012) Primer3, new capabilities and interfaces. *Nucleic Acids Res* 40:e115.
- SantaLucia J (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95:1460–1465.
