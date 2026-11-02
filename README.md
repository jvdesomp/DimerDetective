# DimerDetective

DimerDetective predicts primer dimer formation from the primer3 3′ end stability (`calc_end_stability`) of every primer combination. It comes as a web app that runs entirely in the browser and as a Python command-line script. Both give identical results.

Web app: https://jvdesomp.github.io/DimerDetective/

## What it computes

For each primer pair, DimerDetective computes the free energy (ΔG, kcal/mol) of the 3′ end of each primer annealed to the other primer: `end_left-right_dg` and `end_right-left_dg`. The lower of the two, `min_end_hetero`, is the dimer risk metric used for classification. A more negative ΔG means a more stable dimer and a higher risk. "Left" and "right" refer to the order of the two primers, not to their role.

- **Pair mode** analyzes defined primer pairs.
- **Multiplex mode** tests every primer in a set against every other primer: N primers give N × (N − 1) / 2 heterodimers.

The reaction conditions match the calibration and are fixed, so that results stay comparable across studies: 50 mM Na⁺, 1.5 mM Mg²⁺, 0.6 mM dNTP, 50 nM primer, 37 °C and a maximum loop of 30 nt.

The output is a results table (.xlsx, .tsv or .csv) and an alignment report (.txt). Each alignment diagram is the structure primer3 used to compute the reported ΔG, including bulges and internal loops.

## Training data and performance

The ΔG thresholds were calibrated on 200 RT-qPCR assays: 100 that form primer dimers and 100 clean assays (data courtesy of Biogazelle). The reactions used 250 nM of each primer, SsoAdvanced Universal SYBR Green Supermix, 45 cycles and 60 °C annealing. A positive signal in the no-template control (NTC) defined a primer dimer.

`min_end_hetero` separates the two groups with an area under the ROC curve (AUC) of 0.939.

| `min_end_hetero` threshold | criterion | sensitivity | specificity |
|---|---|---|---|
| < −2.44 kcal/mol | 100 % specificity | 82 % | 100 % |
| < −1.00 kcal/mol | 95 % sensitivity (default) | 95 % | 55 % |
| < −0.27 kcal/mol | 100 % sensitivity | 100 % | 28 % |

| `min_end_hetero` range | risk tier | dimer rate in training data |
|---|---|---|
| < −2.44 kcal/mol | high risk | 100 % (82 dimers, 0 clean) |
| −2.44 to −0.27 kcal/mol | medium risk | 20 % (18 dimers, 72 clean) |
| > −0.27 kcal/mol | low risk | 0 % (0 dimers, 28 clean) |

## Homodimers: not validated

Classification uses `min_end_hetero` only. Homodimer values (`end_left-left_dg`, `end_right-right_dg` and `min_end_homo` in pair mode; `end_self_dg` and the homodimer rows in multiplex mode) are reported for completeness, but they are not validated. In multiplex mode, each primer's self-dimer is reported once, in its homodimer row (`end_self_dg`); heterodimer rows do not repeat it. In the Biogazelle training data, the NTC signal cannot tell homodimers from heterodimers, and predicted homodimers had no predictive value:

- `min_end_homo` gave an AUC of 0.60, against 0.94 for `min_end_hetero`, and added nothing to `min_end_hetero` in a logistic model (p = 0.30).
- 7 of the 100 clean pairs had a strong predicted homodimer (below −2.44 kcal/mol, down to −6.46 kcal/mol) without any NTC signal.

Interpret homodimer values with caution.

## Input

Excel (.xlsx), .csv, .tsv or pasted text. Column names and a header row are optional: sequence columns are recognized by their content, each ID as the column next to its sequence, and an outcome column as a column holding only `dimer` and `clean`; other columns are ignored. Recognized column names (for example `forward_primer`, `left_primer_name`, `primer_seq`, `outcome`) take precedence. The optional outcome column holds the known result of each pair; it is not used in the calculation, but the alignment report compares it with the prediction and the results table repeats it.

## How it works

The thermodynamics are not reimplemented. The C source of primer3's `thal` algorithm (`thal.c`) and its parameter files are taken unchanged from [primer3-py](https://github.com/libnano/primer3-py) 2.4.0 and compiled to WebAssembly with Emscripten:

- `wasm/dd_thal.c` makes the same `thal()` call as `primer3.calc_end_stability()`: alignment type `END1`, mode `THL_GENERAL` and the reaction conditions above. It also asks `thal()` for the secondary structure, so each alignment diagram is the structure behind the reported ΔG.
- `web/dd-core.js` ports the input parsing, scoring metrics and output formats of `python/make_primer_alignments.py`.
- `web/app.js` is the interface. Sequences never leave the browser.

## Python script

Use primer3-py 2.4.0, the release the web app is built from; with it, results are identical to the web app. ΔG comes from `primer3.calc_end_stability`. For the alignment diagrams, the script calls the C function behind it, `thal()`, in the same installed primer3-py, and asks it for the structure as well, so each diagram is the structure behind the reported ΔG.

The script runs on macOS and Linux. On Windows, Python cannot reach `thal()` inside primer3-py and the script stops with a message; use the web app, or run the script in WSL (Windows Subsystem for Linux).

```bash
pip install primer3-py==2.4.0 openpyxl pandas
python python/make_primer_alignments.py --input pairs.xlsx --output results.tsv --threshold -1.00 --type both
python python/make_primer_alignments.py --mode multiplex --input primers.xlsx --output results.xlsx --sort asc
```

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

## How to cite

DimerDetective version 1.2, Jo Vandesompele, Ghent University, https://jvdesomp.github.io/DimerDetective/

Every output file of the web app and the script states the version and build, for example `1.2 (build 20261102-0058)`: in the header of the alignment report and in a final `version` column of the tables. The build is the date and time (Belgian time, YYYYMMDD-hhmm) of the last change; small changes keep version 1.2 and get a new build. The web app shows the build in its footer; the script prints it, also with `--version`.

## Software versions

The results reported here were obtained with:

- DimerDetective 1.2 (build 20261102-0157): web app and Python script
- primer3-py 2.4.0, which includes libprimer3 release 2.6.1: `calc_end_stability` and `thal.c` (the latter compiled unchanged for the web app)
- Emscripten 6.0.13: compilation of `thal.c` to WebAssembly
- Python 3.13.5 with pandas 3.1.0 and openpyxl 3.1.5: Python script (the automated checks also use Python 3.12)
- SheetJS 0.20.3: Excel import and export in the web app
- matplotlib 3.11.2: calibration figures

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

## License and references

GPL-2.0, as primer3 and primer3-py. © Ghent University 2026.

- Untergasser A, Cutcutache I, Koressaar T, Ye J, Faircloth BC, Remm M, Rozen SG (2012) Primer3, new capabilities and interfaces. *Nucleic Acids Res* 40:e115.
- SantaLucia J (1998) A unified view of polymer, dumbbell, and oligonucleotide DNA nearest-neighbor thermodynamics. *PNAS* 95:1460–1465.

## Usage statistics

The web app counts anonymous usage with [GoatCounter](https://www.goatcounter.com) (site code `dimerdetective`): page views, runs per mode and interaction type, the number of pairs or primers per run, and downloads. GoatCounter uses no cookies and stores no personal data. Sequences, primer names and file names are never sent. Nothing is counted when the page is opened from `localhost` or as a local file.

## Contact form

The info tab has a contact form (email address, subject, message; no file upload). Messages are relayed by [Web3Forms](https://web3forms.com) to the developer by email, with the sender's address as reply-to. They contain what the sender enters, plus the sender's IP address, which Web3Forms adds. A hidden field stops simple spam bots. A (quick) reply cannot be guaranteed.
