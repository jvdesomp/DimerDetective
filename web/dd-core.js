/*
 * dd-core.js - DimerDetective computation and report logic.
 *
 * Port of python/make_primer_alignments.py. The thermodynamics run in
 * primer3's own C code (thal.c from primer3-py 2.4.0) compiled to
 * WebAssembly (thal.js); this file only parses input, calls that code,
 * parses the structure and formats the output exactly like the script.
 *
 * Works in the browser (global `DD`) and in Node (module.exports).
 * Copyright (C) 2026 Ghent University. GPL-2.0.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.DD = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // Version (cite this) and build stamp (date and time of the change in
  // Belgian time, YYYYMMDD-hhmm). Must equal VERSION and BUILD in
  // python/make_primer_alignments.py; both go into every output file.
  const VERSION = '1.2';
  const BUILD = '20261102-0951';
  const VERSION_LABEL = `${VERSION} (build ${BUILD})`;

  // Reaction conditions used for calibration (same as PRIMER3_PARAMS).
  const PARAMS = { mv_conc: 50.0, dv_conc: 1.5, dntp_conc: 0.6, dna_conc: 50.0, temp_c: 37.0, max_loop: 30 };

  // Calibrated thresholds (kcal/mol) and risk tiers, materials and methods v1.2.
  const THRESHOLDS = { spec100: -2.44, sens95: -1.00, sens100: -0.27 };
  const DEFAULT_THRESHOLD = THRESHOLDS.sens95;

  // ── WebAssembly engine ──────────────────────────────────────────────────

  let endFn = null, structFn = null, msgFn = null;
  const cache = new Map();

  async function init(createThal) {
    const M = await createThal();
    const rc = M.ccall('dd_init', 'number', ['string'], ['/primer3_config/']);
    if (rc !== 0) throw new Error('could not load primer3 thermodynamic parameters (' + rc + ')');
    endFn = M.cwrap('dd_end_stability', 'number',
      ['string', 'string', 'number', 'number', 'number', 'number', 'number', 'number']);
    structFn = M.cwrap('dd_last_structure', 'string', []);
    msgFn = M.cwrap('dd_last_message', 'string', []);
    cache.clear();
  }

  // 3' end stability of seq1 against seq2: { dg (kcal/mol), structure }.
  // Same call as primer3.calc_end_stability(seq1, seq2, **PRIMER3_PARAMS).
  function endStability(seq1, seq2) {
    const key = seq1 + '|' + seq2;
    let r = cache.get(key);
    if (r) return r;
    const p = PARAMS;
    const dg = endFn(seq1, seq2, p.mv_conc, p.dv_conc, p.dntp_conc, p.dna_conc, p.temp_c, p.max_loop);
    r = msgFn()
      ? { dg: 0.0, structure: '' }          // primer3 error: script falls back to 0.0
      : { dg: dg / 1000, structure: structFn() };
    cache.set(key, r);
    return r;
  }

  function clearCache() { cache.clear(); }

  // ── Structure parsing (port of _parse_ntthal_structure) ─────────────────

  function parseStructure(output) {
    const seqLines = [], strLines = [];
    output.split(/\r?\n/).forEach(line => {
      const s = line.trim();
      const content = () => (line.includes('\t') ? line.slice(line.indexOf('\t') + 1) : line.slice(4));
      if (s.startsWith('SEQ')) seqLines.push(content());
      else if (s.startsWith('STR')) strLines.push(content());
    });
    if (seqLines.length !== 2 || strLines.length !== 2) return null;

    let lines = seqLines.concat(strLines);
    const W = Math.max(...lines.map(l => l.length));
    lines = lines.map(l => l.padEnd(W));

    const paired = [];
    for (let c = 0; c < W; c++) if (lines[1][c] !== ' ' || lines[2][c] !== ' ') paired.push(c);
    const pMin = paired.length ? Math.min(...paired) : W;
    const pMax = paired.length ? Math.max(...paired) : -1;

    const top = [], bot = [];
    for (let c = 0; c < W; c++) {
      const t0 = lines[0][c], t1 = lines[1][c], b2 = lines[2][c], b3 = lines[3][c];
      let t, b;
      if (t1 !== ' ') t = t1;
      else if (t0 === '-') t = (pMin < c && c < pMax) ? '.' : ' ';
      else if (t0 !== ' ') t = t0;
      else t = ' ';
      if (b2 !== ' ') b = b2;
      else if (b3 === '-') b = (pMin < c && c < pMax) ? '.' : ' ';
      else if (b3 !== ' ') b = b3;
      else b = ' ';
      top.push(t); bot.push(b);
    }

    const WC = new Set(['AT', 'TA', 'GC', 'CG']);
    const match = [];
    for (let c = 0; c < W; c++) {
      if (lines[1][c] !== ' ' || lines[2][c] !== ' ') {
        match.push(WC.has(top[c].toUpperCase() + bot[c].toUpperCase()) ? '|' : '.');
      } else {
        match.push(' ');
      }
    }

    const topS = top.join('').trimEnd();
    const matchS = match.join('').trimEnd();
    const botS = bot.join('').trimEnd();

    let consec3p = 0;
    for (let i = matchS.length - 1; i >= 0 && matchS[i] === '|'; i--) consec3p++;

    return {
      top_line: "5'-" + topS + "-3'",
      match_line: '   ' + matchS,
      bot_line: "3'-" + botS + "-5'",
      wc_pairs: (matchS.match(/\|/g) || []).length,
      consec_3p: consec3p,
      method: 'ntthal',
    };
  }

  // ── Per-pair computation (port of compute_pair) ─────────────────────────

  function computePair(seqA, seqB) {
    seqA = seqA.toUpperCase();
    seqB = seqB.toUpperCase();

    const lr = endStability(seqA, seqB);
    const rl = endStability(seqB, seqA);
    const ll = endStability(seqA, seqA);
    const rr = endStability(seqB, seqB);

    const which = lr.dg <= rl.dg ? 'A' : 'B';
    const dom = which === 'A' ? lr : rl;
    const alt = which === 'A' ? rl : lr;

    return {
      'end_left-right_dg': lr.dg,
      'end_right-left_dg': rl.dg,
      'min_end_hetero': Math.min(lr.dg, rl.dg),
      'end_left-left_dg': ll.dg,
      'end_right-right_dg': rr.dg,
      'min_end_homo': Math.min(ll.dg, rr.dg),
      which: which,
      aln: parseStructure(dom.structure),
      aln_alt: seqA === seqB ? null : parseStructure(alt.structure),
    };
  }

  // ── Number formatting identical to Python ───────────────────────────────

  // f"{x:{width}.{prec}f}"
  function fmt(x, prec, width) {
    let s = x.toFixed(prec);
    if (Object.is(x, -0) || (x < 0 && /^0\.?0*$/.test(s))) s = '-' + s;
    return width ? s.padStart(width) : s;
  }

  // str(round(x, 4)) as written by pandas
  function pyRound4(x) {
    const v = Number(x.toFixed(4));
    if (Object.is(v, -0) || (v === 0 && x < 0)) return '-0.0';
    return Number.isInteger(v) ? v.toFixed(1) : String(v);
  }

  // ── Text output (ports of _aln_block and the format_* functions) ────────

  const SEP = '─'.repeat(72);

  function alnBlock(aln, topId, botId, indent) {
    indent = indent === undefined ? '  ' : indent;
    if (!aln) return [`${indent}alignment of 3' ${topId} vs ${botId}: [no significant structure]`];
    return [
      `${indent}alignment of 3' ${topId} vs ${botId} (${aln.wc_pairs} WC pairs, ${aln.consec_3p} consecutive at 3' end)`,
      `${indent}${aln.top_line}`,
      `${indent}${aln.match_line}`,
      `${indent}${aln.bot_line}`,
    ];
  }

  function directionIds(result, idA, idB) {
    return result.which === 'A'
      ? { domTop: idA, domBot: idB, altTop: idB, altBot: idA, dominant: idA }
      : { domTop: idB, domBot: idA, altTop: idA, altBot: idB, dominant: idB };
  }

  function formatPair(row, result, threshold) {
    const f = row.fwd_name, r = row.rev_name;
    const lines = [SEP, `${f} / ${r}  [heterodimer]`, ''];
    const call = result.min_end_hetero < threshold ? 'dimer' : 'clean';
    if (row.outcome !== undefined && row.outcome !== null && row.outcome !== '') {
      const label = String(row.outcome).toLowerCase();
      const correct = call === label ? '✓' : '✗ (label mismatch)';
      lines.push(`  predicted: ${call}   label: ${label}   ${correct}`);
    } else {
      lines.push(`  predicted: ${call}`);
    }
    const d = directionIds(result, f, r);
    lines.push(`  end_left-right_dg = ${fmt(result['end_left-right_dg'], 3, 7)} kcal/mol  (3' of ${f} vs ${r})`);
    lines.push(`  end_right-left_dg = ${fmt(result['end_right-left_dg'], 3, 7)} kcal/mol  (3' of ${r} vs ${f})`);
    lines.push(`  min_end_hetero = ${fmt(result.min_end_hetero, 3, 7)} kcal/mol  ← dominant (${d.dominant} 3'-end)`);
    lines.push('');
    lines.push(`  ${f}: 5'-${row.fwd}-3'`);
    lines.push(`  ${r}: 5'-${row.rev}-3'`);
    lines.push('');
    lines.push(...alnBlock(result.aln, d.domTop, d.domBot));
    lines.push('');
    lines.push(...alnBlock(result.aln_alt, d.altTop, d.altBot));
    lines.push('');
    return lines;
  }

  function formatHomodimerEntry(pid, seq, result) {
    return [
      SEP, `${pid} / ${pid}  [homodimer]`, '',
      `  min_end_homo = ${fmt(result.min_end_hetero, 3, 7)} kcal/mol`, '',
      `  ${pid}: 5'-${seq}-3'`, '',
      ...alnBlock(result.aln, pid, pid), '',
    ];
  }

  function formatMultiplexEntry(e, threshold) {
    const { type, idA, seqA, idB, seqB, result } = e;
    const call = result.min_end_hetero < threshold ? 'dimer' : 'clean';
    const lines = [SEP];
    if (type === 'homodimer') {
      lines.push(`${idA} / ${idA}  [homodimer]`, '');
      lines.push(`  predicted: ${call}`);
      lines.push(`  end_self_dg = ${fmt(result['end_left-left_dg'], 3, 7)} kcal/mol`, '');
      lines.push(`  ${idA}: 5'-${seqA}-3'`, '');
      lines.push(...alnBlock(result.aln, idA, idA));
    } else {
      const d = directionIds(result, idA, idB);
      lines.push(`${idA} / ${idB}  [heterodimer]`, '');
      lines.push(`  predicted: ${call}`);
      lines.push(`  end_left-right_dg = ${fmt(result['end_left-right_dg'], 3, 7)} kcal/mol  (3' of ${idA} vs ${idB})`);
      lines.push(`  end_right-left_dg = ${fmt(result['end_right-left_dg'], 3, 7)} kcal/mol  (3' of ${idB} vs ${idA})`);
      lines.push(`  min_end_hetero    = ${fmt(result.min_end_hetero, 3, 7)} kcal/mol  ← dominant (${d.dominant} 3'-end)`);
      lines.push('');
      lines.push(`  ${idA}: 5'-${seqA}-3'`);
      lines.push(`  ${idB}: 5'-${seqB}-3'`, '');
      lines.push(...alnBlock(result.aln, d.domTop, d.domBot), '');
      lines.push(...alnBlock(result.aln_alt, d.altTop, d.altBot));
    }
    lines.push('');
    return lines;
  }

  const SORT_LABELS = {
    asc: 'ascending (most stable dimer first)',
    desc: 'descending (least stable dimer first)',
    input: 'input order',
  };

  // ── Input parsing (port of find_col / _read_primer_list / run_pairs) ────

  const DNA_PAT = /^[ATGCUNatgcun-]+$/;
  const ACGT = /^[ACGT]+$/;

  function looksLikeDna(val, minLen) {
    const s = String(val === undefined || val === null ? '' : val).trim();
    return s.length >= (minLen || 6) && DNA_PAT.test(s);
  }

  function isMissing(v) {
    if (v === undefined || v === null) return true;
    const s = String(v).trim();
    return s === '' || ['nan', 'none'].includes(s.toLowerCase());
  }

  function dnaFraction(rows, c) {
    const vals = rows.map(r => r[c]).filter(v => !isMissing(v)).map(v => String(v).trim());
    return vals.length ? vals.filter(v => looksLikeDna(v)).length / vals.length : 0;
  }

  function findCol(columns, candidates) {
    const lower = new Map();
    columns.forEach((c, i) => { if (!lower.has(String(c).toLowerCase())) lower.set(String(c).toLowerCase(), i); });
    for (const cand of candidates) {
      const k = cand.toLowerCase();
      if (lower.has(k)) return lower.get(k);
    }
    return null;
  }

  // Split one line on `sep`, honoring "quoted" cells.
  function splitQuoted(line, sep) {
    const cells = [];
    let cur = '', inQ = false;
    for (let i = 0; i < line.length; i++) {
      const ch = line[i];
      if (inQ) {
        if (ch === '"' && line[i + 1] === '"') { cur += '"'; i++; }
        else if (ch === '"') inQ = false;
        else cur += ch;
      } else if (ch === '"' && cur.trim() === '') { inQ = true; cur = ''; }
      else if (ch === sep) { cells.push(cur); cur = ''; }
      else cur += ch;
    }
    cells.push(cur);
    return cells;
  }

  // Split pasted text into rows of cells (tab, comma, semicolon or spaces).
  function textToRows(text) {
    const rows = [];
    text.split(/\r?\n/).forEach(line => {
      if (!line.trim() || line.trim().startsWith('#')) return;
      const unquoted = line.replace(/"[^"]*"/g, '');
      let cells;
      if (unquoted.includes('\t')) cells = splitQuoted(line, '\t');
      else if (unquoted.includes(',')) cells = splitQuoted(line, ',');
      else if (unquoted.includes(';')) cells = splitQuoted(line, ';');
      else cells = line.trim().split(/\s+/);
      rows.push(cells.map(c => c.trim()));
    });
    return rows;
  }

  // A header row is present when the first row holds no sequence.
  function splitHeader(rows) {
    if (!rows.length) return { columns: [], data: [] };
    const width = Math.max(...rows.map(r => r.length));
    const pad = r => { const x = r.slice(); while (x.length < width) x.push(''); return x; };
    const first = pad(rows[0]);
    const hasHeader = !first.some(c => looksLikeDna(c));
    const columns = hasHeader ? first.map((c, i) => c || 'column_' + (i + 1)) : first.map((_, i) => i);
    const data = (hasHeader ? rows.slice(1) : rows).map(pad);
    return { columns, data, hasHeader };
  }

  // Column names are optional. Sequence columns are recognized by their
  // content, each ID column as the non-sequence column next to its sequence
  // column, and the outcome column as a column holding only "dimer" and
  // "clean". Recognized names take precedence. Mirrors _pair_columns,
  // _read_pairs and _read_primer_list in the script; both must stay identical.
  const PAIR_COLS = {
    fwdName: ['forward_primer_name', 'fwd_name', 'fwd_id', 'forward_name', 'forward_id',
      'left_primer_name', 'left_name', 'left_id', 'left-id', 'Forward Primer Name'],
    fwdSeq: ['forward_primer', 'forward_primer_seq', 'forward_primer_sequence',
      'forward_sequence', 'fwd', 'fwd_seq', 'left_primer', 'left_primer_seq', 'left_primer_sequence',
      'left_sequence', 'left-sequence', 'left-primer', 'Forward Primer Seq', 'Forward Primer'],
    revName: ['reverse_primer_name', 'rev_name', 'rev_id', 'reverse_name', 'reverse_id',
      'right_primer_name', 'right_name', 'right_id', 'right-id', 'Reverse Primer Name'],
    revSeq: ['reverse_primer', 'reverse_primer_seq', 'reverse_primer_sequence',
      'reverse_sequence', 'rev', 'rev_seq', 'right_primer', 'right_primer_seq', 'right_primer_sequence',
      'right_sequence', 'right-sequence', 'right-primer', 'Reverse Primer Seq', 'Reverse Primer'],
    outcome: ['outcome', 'label', 'class', 'result', 'status'],
  };
  const LIST_COLS = {
    seq: ['primer_seq', 'primer_sequence', 'sequence', 'seq'],
    id: ['primer_name', 'primer_id', 'name', 'id'],
  };
  const OUTCOME_VALUES = new Set(['dimer', 'clean']);

  function cleanSeq(v) { return String(v === undefined || v === null ? '' : v).replace(/\s+/g, ''); }

  function cellStr(v) {
    if (typeof v === 'number' && Number.isInteger(v)) return String(v);
    return String(v === undefined || v === null ? '' : v).trim();
  }

  // Identifier text of a cell ('1.0' -> '1'), or null when empty.
  function cellId(v) {
    if (isMissing(v)) return null;
    let s = cellStr(v);
    if (/^-?\d+\.0+$/.test(s)) s = String(parseInt(s, 10));
    return s;
  }

  function dnaColumns(columns, data) {
    return columns.map((_, i) => i).filter(i => dnaFraction(data, i) >= 0.8);
  }

  function outcomeColumn(columns, data, exclude) {
    const named = findCol(columns, PAIR_COLS.outcome);
    if (named !== null) return named;
    for (let c = 0; c < columns.length; c++) {
      if (exclude.has(c)) continue;
      const vals = data.map(r => r[c]).filter(v => !isMissing(v)).map(v => String(v).trim().toLowerCase());
      if (vals.length && vals.every(v => OUTCOME_VALUES.has(v))) return c;
    }
    return null;
  }

  // The non-sequence column directly left of seqCol, else directly right.
  function idColumnNextTo(columns, data, seqCol, taken) {
    for (const j of [seqCol - 1, seqCol + 1]) {
      if (j >= 0 && j < columns.length && !taken.has(j) && dnaFraction(data, j) < 0.5) return j;
    }
    return null;
  }

  function pairColumns(columns, data) {
    let fName = findCol(columns, PAIR_COLS.fwdName);
    let fSeq = findCol(columns, PAIR_COLS.fwdSeq);
    let rName = findCol(columns, PAIR_COLS.revName);
    let rSeq = findCol(columns, PAIR_COLS.revSeq);
    if (fSeq === null || rSeq === null) {
      const dna = dnaColumns(columns, data).filter(c => c !== fSeq && c !== rSeq);
      if (fSeq === null && dna.length) fSeq = dna.shift();
      if (rSeq === null && dna.length) rSeq = dna.shift();
    }
    if (fSeq === null || rSeq === null) return null;
    const out = outcomeColumn(columns, data, new Set([fSeq, rSeq, fName, rName].filter(c => c !== null)));
    const taken = new Set([fSeq, rSeq, out, fName, rName].filter(c => c !== null));
    if (fName === null) { fName = idColumnNextTo(columns, data, fSeq, taken); if (fName !== null) taken.add(fName); }
    if (rName === null) rName = idColumnNextTo(columns, data, rSeq, taken);
    return { fName, fSeq, rName, rSeq, out };
  }

  // Pair mode: returns { rows: [{fwd_name, fwd, rev_name, rev, outcome}], errors, hasOutcome }
  function parsePairs(rawRows) {
    const { columns, data } = splitHeader(rawRows);
    const errors = [];
    const cols = pairColumns(columns, data);
    if (!cols) {
      errors.push('could not find two primer sequence columns (left and right)');
      return { rows: [], errors, hasOutcome: false };
    }
    const { fName, fSeq, rName, rSeq, out } = cols;
    const rows = [];
    data.forEach((r, i) => {
      if (isMissing(r[fSeq]) && isMissing(r[rSeq])) return;
      const fwd = cellStr(r[fSeq]), rev = cellStr(r[rSeq]);
      if (!ACGT.test(cleanSeq(fwd).toUpperCase()) || !ACGT.test(cleanSeq(rev).toUpperCase())) {
        errors.push(`row ${i + 1}: sequences may only contain A, C, G and T`);
        return;
      }
      rows.push({
        fwd_name: (fName !== null ? cellId(r[fName]) : null) || `FWD_${i + 1}`,
        rev_name: (rName !== null ? cellId(r[rName]) : null) || `REV_${i + 1}`,
        fwd: cleanSeq(fwd),
        rev: cleanSeq(rev),
        outcome: out !== null && !isMissing(r[out]) ? cellStr(r[out]) : null,
      });
    });
    return { rows, errors, hasOutcome: out !== null };
  }

  // Multiplex mode: returns { primers: [[id, seq]], errors, notes }
  function parsePrimerList(rawRows) {
    const { columns, data } = splitHeader(rawRows);
    const errors = [];
    const nDigits = String(data.length).length;
    const autoId = n => 'P' + String(n).padStart(nDigits, '0');

    let seqCol = findCol(columns, LIST_COLS.seq);
    if (seqCol === null) {
      const dna = dnaColumns(columns, data);
      // A table of primer pairs (two sequence columns): use all primers of both.
      if (dna.length >= 2 || (findCol(columns, PAIR_COLS.fwdSeq) !== null && findCol(columns, PAIR_COLS.revSeq) !== null)) {
        const pairs = parsePairs(rawRows);
        const found = [];
        pairs.rows.forEach(p => { found.push([p.fwd_name, p.fwd.toUpperCase()], [p.rev_name, p.rev.toUpperCase()]); });
        const { primers, notes } = uniquePrimers(found);
        notes.unshift(`primer pair table: primers from both sequence columns are used (${found.length} primers)`);
        return { primers, errors: pairs.errors, notes };
      }
      if (!dna.length) {
        errors.push('could not identify a primer sequence column; add a header named primer_seq');
        return { primers: [], errors, notes: [] };
      }
      seqCol = dna[0];
    }
    let idCol = findCol(columns, LIST_COLS.id);
    if (idCol === null) idCol = idColumnNextTo(columns, data, seqCol, new Set([seqCol]));

    const primers = [];
    data.forEach((r, i) => {
      if (isMissing(r[seqCol])) return;
      const seq = cleanSeq(r[seqCol]).toUpperCase();
      if (!ACGT.test(seq)) { errors.push(`row ${i + 1}: sequences may only contain A, C, G and T`); return; }
      const pid = (idCol !== null ? cellId(r[idCol]) : null) || autoId(primers.length + 1);
      primers.push([pid, seq]);
    });
    return { ...uniquePrimers(primers), errors };
  }

  // Keep the first primer of each sequence; report merged duplicates
  // (same wording as _unique_primers in the script).
  function uniquePrimers(list) {
    const kept = [], seen = new Map(), dups = [];
    list.forEach(([pid, seq]) => {
      if (seen.has(seq)) dups.push([pid, seen.get(seq)]);
      else { seen.set(seq, pid); kept.push([pid, seq]); }
    });
    const notes = [];
    if (dups.length) {
      notes.push(`${list.length} primers in the input, ${kept.length} unique sequences; ` +
        `duplicates merged (first ID kept): ${dups.map(([d, k]) => `${d} = ${k}`).join(', ')}`);
    }
    return { primers: kept, notes };
  }

  // ── Analyses ────────────────────────────────────────────────────────────

  // Each analysis returns { mode, type, entries, ... }; entries hold
  // { type, idA, seqA, idB, seqB, result, n, outcome }.

  function pairEntries(rows, type, onProgress) {
    const entries = [];
    const pairResults = rows.map((row, i) => {
      const result = computePair(row.fwd, row.rev);
      if (onProgress) onProgress(i + 1, rows.length);
      return result;
    });
    if (type === 'heterodimers' || type === 'both') {
      rows.forEach((row, i) => entries.push({
        type: 'heterodimer', n: i + 1, idA: row.fwd_name, seqA: row.fwd, idB: row.rev_name, seqB: row.rev,
        outcome: row.outcome, result: pairResults[i],
      }));
    }
    if (type === 'homodimers' || type === 'both') {
      const seen = new Set();
      rows.forEach(row => {
        [[row.fwd_name, row.fwd], [row.rev_name, row.rev]].forEach(([pid, seq]) => {
          if (seen.has(pid)) return;
          seen.add(pid);
          entries.push({ type: 'homodimer', n: entries.length + 1, idA: pid, seqA: seq, idB: pid, seqB: seq, result: computePair(seq, seq) });
        });
      });
    }
    return { pairResults, entries };
  }

  function multiplexJobs(primers, type) {
    const jobs = [];
    if (type === 'heterodimers' || type === 'both') {
      for (let i = 0; i < primers.length; i++)
        for (let j = i + 1; j < primers.length; j++)
          jobs.push(['heterodimer', primers[i], primers[j]]);
    }
    if (type === 'homodimers' || type === 'both') primers.forEach(p => jobs.push(['homodimer', p, p]));
    return jobs;
  }

  function runJob(job, n) {
    const [type, [idA, seqA], [idB, seqB]] = job;
    return { type, n, idA, seqA, idB, seqB, result: computePair(seqA, seqB) };
  }

  function sortEntries(entries, sort) {
    const e = entries.slice();
    if (sort === 'asc') e.sort((a, b) => a.result.min_end_hetero - b.result.min_end_hetero);
    else if (sort === 'desc') e.sort((a, b) => b.result.min_end_hetero - a.result.min_end_hetero);
    return e;
  }

  // ── Reports (exact script output) ───────────────────────────────────────

  function pairsReport(rows, pairResults, threshold, type) {
    const flagged = pairResults.filter(r => r.min_end_hetero < threshold).length;
    const lines = [
      '='.repeat(72), 'DimerDetective alignment output', '',
      `version   : ${VERSION_LABEL}`,
      'mode      : pairs',
      `type      : ${type}`,
      `threshold : ${fmt(threshold, 2)} kcal/mol → dimer`,
      `pairs     : ${rows.length}   |   flagged: ${flagged}`,
      '| = Watson-Crick pair', '. = bulge on opposite strand',
      '='.repeat(72), '',
    ];
    if (type === 'heterodimers' || type === 'both') {
      rows.forEach((row, i) => lines.push(...formatPair(row, pairResults[i], threshold)));
    }
    if (type === 'homodimers' || type === 'both') {
      lines.push('', '='.repeat(72), 'homodimers', '='.repeat(72), '');
      const seen = new Set();
      rows.forEach(row => {
        [[row.fwd_name, row.fwd], [row.rev_name, row.rev]].forEach(([pid, seq]) => {
          if (seen.has(pid)) return;
          seen.add(pid);
          lines.push(...formatHomodimerEntry(pid, seq, computePair(seq, seq)));
        });
      });
    }
    return lines.join('\n');
  }

  function pairsTable(rows, pairResults, threshold, hasOutcome) {
    return rows.map((row, i) => {
      const r = pairResults[i];
      const flagged = r.min_end_hetero < threshold;
      const t = {
        fwd_name: row.fwd_name, fwd_seq: row.fwd, rev_name: row.rev_name, rev_seq: row.rev,
        'end_left-right_dg': r['end_left-right_dg'], 'end_right-left_dg': r['end_right-left_dg'],
        min_end_hetero: r.min_end_hetero,
        'end_left-left_dg': r['end_left-left_dg'], 'end_right-right_dg': r['end_right-right_dg'],
        min_end_homo: r.min_end_homo,
        flagged: flagged ? 'YES' : 'NO', predicted: flagged ? 'dimer' : 'clean',
      };
      if (hasOutcome) t.outcome = row.outcome === null ? '' : row.outcome;
      t.version = VERSION_LABEL;
      return t;
    });
  }

  function multiplexReport(nPrimers, sorted, threshold, type, sort, notes) {
    const flagged = sorted.filter(e => e.result.min_end_hetero < threshold).length;
    const lines = [
      '='.repeat(72), 'DimerDetective alignment output', '',
      `version   : ${VERSION_LABEL}`,
      'mode      : multiplex',
      `threshold : ${fmt(threshold, 2)} kcal/mol → dimer`,
      `primers   : ${nPrimers}   |   pairs: ${sorted.length}   |   flagged: ${flagged}`,
      ...(notes || []).map(n => `note      : ${n}`),
      `type      : ${type}`,
      `sorted by : ${SORT_LABELS[sort] || sort}`,
      '| = Watson-Crick pair', '. = bulge on opposite strand',
      '='.repeat(72), '',
    ];
    sorted.forEach(e => lines.push(...formatMultiplexEntry(e, threshold)));
    return lines.join('\n');
  }

  function multiplexTable(sorted, threshold) {
    return sorted.map(e => {
      const r = e.result;
      const homo = e.type === 'homodimer';
      return {
        type: e.type, id_A: e.idA, seq_A: e.seqA, id_B: e.idB, seq_B: e.seqB,
        // heterodimer rows: the three heterodimer values; homodimer rows:
        // one self-dimer value (end_self_dg), the others left empty
        'end_left-right_dg': homo ? null : r['end_left-right_dg'],
        'end_right-left_dg': homo ? null : r['end_right-left_dg'],
        min_end_hetero: homo ? null : r.min_end_hetero,
        end_self_dg: homo ? r['end_left-left_dg'] : null,
        flagged: r.min_end_hetero < threshold ? 'YES' : 'NO',
        version: VERSION_LABEL,
      };
    });
  }

  const DG_COLS = ['end_left-right_dg', 'end_right-left_dg', 'min_end_hetero', 'end_left-left_dg', 'end_right-right_dg', 'min_end_homo', 'end_self_dg'];

  // Delimited text like pandas to_csv (values rounded to 4 decimals).
  function tableToText(table, sep) {
    if (!table.length) return '';
    const cols = Object.keys(table[0]);
    const cell = v => {
      let s = String(v);
      if (sep === ',' && /[",\n\r]/.test(s)) s = '"' + s.replace(/"/g, '""') + '"';
      return s;
    };
    const lines = [cols.join(sep)];
    table.forEach(row => lines.push(cols.map(c => cell(row[c] === null ? '' : DG_COLS.includes(c) ? pyRound4(row[c]) : row[c])).join(sep)));
    return lines.join('\n') + '\n';
  }

  // Rows for spreadsheet export (numbers rounded to 4 decimals).
  function tableToAoa(table) {
    if (!table.length) return [];
    const cols = Object.keys(table[0]);
    return [cols].concat(table.map(row => cols.map(c => (row[c] === null ? null : DG_COLS.includes(c) ? Number(row[c].toFixed(4)) : row[c]))));
  }

  function riskTier(dg) {
    if (dg < THRESHOLDS.spec100) return 'high';
    if (dg <= THRESHOLDS.sens100) return 'medium';
    return 'low';
  }

  return {
    VERSION, BUILD, VERSION_LABEL, PARAMS, THRESHOLDS, DEFAULT_THRESHOLD,
    init, endStability, clearCache, parseStructure, computePair,
    textToRows, parsePairs, parsePrimerList,
    pairEntries, multiplexJobs, runJob, sortEntries,
    pairsReport, pairsTable, multiplexReport, multiplexTable,
    tableToText, tableToAoa, riskTier, fmt, alnBlock, SORT_LABELS,
  };
});
