#!/usr/bin/env node
// Command-line front end for web/dd-core.js with the same arguments as
// python/make_primer_alignments.py. Used to check that the web version
// produces byte-identical output. Input: tab-separated text (.tsv/.txt/.csv).
'use strict';
const fs = require('fs');
const path = require('path');
const DD = require('../web/dd-core.js');
const createThal = require('../web/thal.js');

function parseArgs(argv) {
  const a = { mode: 'pairs', threshold: DD.DEFAULT_THRESHOLD, sort: 'asc', type: 'heterodimers', alignments: null };
  for (let i = 0; i < argv.length; i += 2) {
    const k = argv[i].replace(/^--/, ''), v = argv[i + 1];
    a[k] = k === 'threshold' ? parseFloat(v) : v;
  }
  if (!a.input || !a.output) throw new Error('usage: dd-cli.js --input FILE --output FILE [--mode pairs|multiplex] ...');
  if (!a.alignments) a.alignments = a.output.replace(/\.[^.]*$/, '') + '_alignments.txt';
  return a;
}

(async () => {
  const args = parseArgs(process.argv.slice(2));
  await DD.init(createThal);
  const rows = DD.textToRows(fs.readFileSync(args.input, 'utf8'));

  if (args.mode === 'pairs') {
    const parsed = DD.parsePairs(rows);
    if (parsed.errors.length) console.error(parsed.errors.join('\n'));
    const { pairResults } = DD.pairEntries(parsed.rows, 'heterodimers');
    const table = DD.pairsTable(parsed.rows, pairResults, args.threshold, parsed.hasOutcome);
    fs.writeFileSync(args.output, DD.tableToText(table, '\t'));
    fs.writeFileSync(args.alignments, DD.pairsReport(parsed.rows, pairResults, args.threshold, args.type));
    console.log(`Written ${parsed.rows.length} pairs → ${path.basename(args.output)}`);
  } else {
    const parsed = DD.parsePrimerList(rows);
    if (parsed.errors.length) console.error(parsed.errors.join('\n'));
    parsed.notes.forEach(n => console.log(`Note: ${n}`));
    const jobs = DD.multiplexJobs(parsed.primers, args.type);
    const entries = jobs.map((j, i) => DD.runJob(j, i + 1));
    const sorted = DD.sortEntries(entries, args.sort);
    fs.writeFileSync(args.output, DD.tableToText(DD.multiplexTable(sorted, args.threshold), '\t'));
    fs.writeFileSync(args.alignments, DD.multiplexReport(parsed.primers.length, sorted, args.threshold, args.type, args.sort, parsed.notes));
    console.log(`Written ${sorted.length} rows → ${path.basename(args.output)}`);
  }
})().catch(e => { console.error(e.message); process.exit(1); });
