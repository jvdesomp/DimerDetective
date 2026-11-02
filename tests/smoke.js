#!/usr/bin/env node
// Quick check of the WebAssembly build against values from
// primer3.calc_end_stability (primer3-py 2.4.0) at the calibration
// conditions. Runs in CI before deployment.
'use strict';
const DD = require('../web/dd-core.js');
const createThal = require('../web/thal.js');

// [seq1, seq2, dG in cal/mol from primer3-py]
const CASES = [
  ['CTCTATGTGTTGGAGTCATGAAG', 'CATTCCTTCTGAGGCTCATT', -2651.420102457284],
  ['CATTCCTTCTGAGGCTCATT', 'CTCTATGTGTTGGAGTCATGAAG', 95.91377107841812],
  ['GGATCTGACTGTTCAGCAA', 'GAGGAACATCTGTATTGTAGGTT', -1671.332457846264],
  ['ACAACATCGTAAGGAGGAATG', 'ACAACATCGTAAGGAGGAATG', -516.1395603004385],
  ['AGATGTGAGTTATAAATGAAGGTT', 'CCAGGAGGAAAGTTTCAGG', 383.9138553920976],
  ['GCTAGCTAGGTCAACG', 'CGTTGACCTAGC', -12445.648518143606],
];

const ALIGNMENT = [
  "5'-CTCTATGTGTTGGAGTCAT.GAAG-3'",
  '               ||| |   ||||',
  "3'-         TTACTCGGAGTCTTCCTTAC-5'",
];

(async () => {
  await DD.init(createThal);
  let failed = 0;
  for (const [a, b, ref] of CASES) {
    const dg = DD.endStability(a, b).dg * 1000;
    const ok = Math.abs(dg - ref) < 1e-3;
    if (!ok) failed++;
    console.log(`${ok ? 'ok  ' : 'FAIL'} ${a} vs ${b}: ${dg.toFixed(4)} (primer3-py ${ref.toFixed(4)})`);
  }
  const aln = DD.computePair(CASES[0][0], CASES[0][1]).aln;
  const lines = [aln.top_line, aln.match_line, aln.bot_line];
  const alnOk = lines.join('\n') === ALIGNMENT.join('\n');
  if (!alnOk) failed++;
  console.log(`${alnOk ? 'ok  ' : 'FAIL'} alignment\n${lines.join('\n')}`);
  // the script and the web app must carry the same version and build
  const py = require('fs').readFileSync(require('path').join(__dirname, '..', 'python', 'make_primer_alignments.py'), 'utf8');
  const pyVersion = (py.match(/^VERSION = '([^']+)'/m) || [])[1];
  const pyBuild = (py.match(/^BUILD = '([^']+)'/m) || [])[1];
  const verOk = pyVersion === DD.VERSION && pyBuild === DD.BUILD && /^\d{8}-\d{4}$/.test(DD.BUILD);
  if (!verOk) failed++;
  console.log(`${verOk ? 'ok  ' : 'FAIL'} version: web ${DD.VERSION_LABEL}, script ${pyVersion} (build ${pyBuild})`);
  if (failed) { console.error(`${failed} check(s) failed`); process.exit(1); }
  console.log('all checks passed');
})().catch(e => { console.error(e); process.exit(1); });
