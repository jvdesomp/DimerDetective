/*
 * app.js - DimerDetective user interface. Computation and report formats
 * live in dd-core.js; this file handles input, progress, display and
 * downloads. Copyright (C) 2026 Ghent University. GPL-2.0.
 */
'use strict';

const MAX_PAIRS = 1000;
const MAX_PRIMERS = 200;
const DISPLAY_LIMIT = 500;      // rows rendered in the table; downloads hold all
const CHUNK = 25;               // computations between screen updates
const GOATCOUNTER = 'dimerdetective';   // GoatCounter site code (https://CODE.goatcounter.com); empty = no counting
const WEB3FORMS_KEY = '4e58b612-17f4-4b78-8b83-2975281ebd71';   // Web3Forms access key (public by design); empty = contact form disabled

const EXAMPLES = {
  pairs: [
    'left_id\tleft_sequence\tright_id\tright_sequence',
    'pair1-F\tCTCTATGTGTTGGAGTCATGAAG\tpair1-R\tCATTCCTTCTGAGGCTCATT',
    'pair2-F\tGGATCTGACTGTTCAGCAA\tpair2-R\tGAGGAACATCTGTATTGTAGGTT',
    'pair3-F\tAGATGTGAGTTATAAATGAAGGTT\tpair3-R\tCCAGGAGGAAAGTTTCAGG',
  ].join('\n'),
  multiplex: [
    'primer_id\tprimer_seq',
    'pair1-F\tCTCTATGTGTTGGAGTCATGAAG',
    'pair1-R\tCATTCCTTCTGAGGCTCATT',
    'pair2-F\tGGATCTGACTGTTCAGCAA',
    'pair2-R\tGAGGAACATCTGTATTGTAGGTT',
    'pair3-F\tAGATGTGAGTTATAAATGAAGGTT',
    'pair3-R\tCCAGGAGGAAAGTTTCAGG',
  ].join('\n'),
};

const $ = id => document.getElementById(id);
const state = {
  mode: 'pairs',          // input tab
  texts: { pairs: '', multiplex: '' },
  stems: { pairs: null, multiplex: null },
  run: null,              // last analysis
  busy: false,
};

// ── Helpers ───────────────────────────────────────────────────────────────

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;').replace(/'/g, '&#39;');
}

function getThreshold() {
  const sel = $('threshold-preset');
  if (sel.value !== 'custom') return parseFloat(sel.value);
  const v = parseFloat(($('threshold-custom').value || '').replace(',', '.').replace('−', '-'));
  return Number.isFinite(v) ? v : DD.DEFAULT_THRESHOLD;
}

const fmt2 = x => (x === null ? '' : DD.fmt(x, 2));

function riskBadge(dg) {
  const tier = DD.riskTier(dg);
  return `<span class="badge badge-${tier}">${tier === 'none' ? 'not computed' : tier + ' risk'}</span>`;
}

function setMessages(lines) {
  $('messages').textContent = lines.length ? lines.slice(0, 12).join('\n') + (lines.length > 12 ? `\n(${lines.length - 12} more)` : '') : '';
}

function download(content, filename, type) {
  const blob = content instanceof Blob ? content : new Blob([content], { type: type || 'text/plain;charset=utf-8' });
  const a = document.createElement('a');
  a.href = URL.createObjectURL(blob);
  a.download = filename;
  document.body.appendChild(a);
  a.click();
  setTimeout(() => { URL.revokeObjectURL(a.href); a.remove(); }, 1000);
}

function stem() {
  return (state.run && state.run.stem) || 'dimerdetective';
}

// Replace the bulge character on the strands (not in the match string).
function withGap(line, gap) {
  if (gap === '.' || !(line.startsWith("5'-") || line.startsWith("3'-"))) return line;
  return line.slice(0, 3) + line.slice(3, -3).replace(/\./g, gap) + line.slice(-3);
}

function reportWithGap(text, gap) {
  if (gap === '.') return text;
  return text.split('\n').map(l => {
    const m = l.match(/^(\s*)(.*)$/);
    return m[1] + withGap(m[2], gap);
  }).join('\n');
}

// ── Input ─────────────────────────────────────────────────────────────────

function switchMode(mode) {
  state.texts[state.mode] = $('input-text').value;
  state.mode = mode;
  $('input-text').value = state.texts[mode];
  ['pairs', 'multiplex'].forEach(m => {
    $('tab-' + m).classList.toggle('active', m === mode);
    $('hint-' + m).classList.toggle('hidden', m !== mode);
  });
  $('run-btn').innerHTML = '&#9654; run ' + (mode === 'pairs' ? 'pair' : 'multiplex') + ' mode';
  $('input-text').placeholder = EXAMPLES[mode];
  setMessages([]);
}

function cellText(v) {
  return String(v === undefined || v === null ? '' : v).replace(/\s+/g, ' ').trim();
}

async function loadFile(file) {
  if (!file) return;
  const ext = file.name.split('.').pop().toLowerCase();
  let text;
  if (ext === 'xlsx' || ext === 'xls') {
    if (typeof XLSX === 'undefined') {
      setMessages(['Excel files need the SheetJS library, which could not be loaded. Save the sheet as .csv or .tsv, or paste the cells.']);
      return;
    }
    const wb = XLSX.read(await file.arrayBuffer(), { type: 'array' });
    const rows = XLSX.utils.sheet_to_json(wb.Sheets[wb.SheetNames[0]], { header: 1, raw: false, defval: '' });
    text = rows.filter(r => r.some(c => cellText(c))).map(r => r.map(cellText).join('\t')).join('\n');
  } else {
    text = await file.text();
  }
  $('input-text').value = text;
  state.texts[state.mode] = text;
  state.stems[state.mode] = file.name.replace(/\.[^.]*$/, '');
  setMessages([`loaded ${file.name}`]);
}

// ── Anonymous usage counts ────────────────────────────────────────────────
// Events carry the mode, interaction type and numbers only, never sequences,
// primer IDs or file names. GoatCounter drops more than 4 requests per second
// from one address, so events go out one by one.

const pendingEvents = [];
let eventTimer = 0;

function track(...paths) {
  if (!GOATCOUNTER) return;
  pendingEvents.push(...paths);
  sendEvents();
}

function sendEvents() {
  const gc = window.goatcounter;
  if (eventTimer || !pendingEvents.length || !(gc && gc.count)) return;
  const path = pendingEvents.shift();
  gc.count({ path, title: path, event: true });
  eventTimer = setTimeout(() => { eventTimer = 0; sendEvents(); }, 350);
}

function setupCounting() {
  if (!GOATCOUNTER || !/^https?:$/.test(location.protocol)) return;
  const s = document.createElement('script');
  s.async = true;
  s.src = 'https://gc.zgo.at/count.js';
  s.dataset.goatcounter = `https://${GOATCOUNTER}.goatcounter.com/count`;
  // count.js sends the page view itself on load; events follow after it
  s.onload = () => { eventTimer = setTimeout(() => { eventTimer = 0; sendEvents(); }, 350); };
  document.head.appendChild(s);
}

// ── Running ───────────────────────────────────────────────────────────────

function setProgress(done, total, what) {
  $('progress-bar').style.width = (total ? done / total * 100 : 0) + '%';
  $('progress-status').textContent = `computing ${what}: ${done} / ${total}`;
}

function runChunked(total, step, what) {
  return new Promise(resolve => {
    let done = 0;
    $('progress-card').classList.remove('hidden');
    const tick = () => {
      const end = Math.min(done + CHUNK, total);
      for (; done < end; done++) step(done);
      setProgress(done, total, what);
      if (done < total) setTimeout(tick, 0);
      else { $('progress-card').classList.add('hidden'); resolve(); }
    };
    setTimeout(tick, 0);
  });
}

async function run() {
  if (state.busy) return;
  const text = $('input-text').value;
  state.texts[state.mode] = text;
  const rows = DD.textToRows(text);
  const type = $('type-mode').value;
  const messages = [];
  DD.clearCache();
  state.busy = true;
  $('run-btn').disabled = true;
  try {
    if (state.mode === 'pairs') {
      const parsed = DD.parsePairs(rows);
      messages.push(...parsed.errors);
      if (!parsed.rows.length) { messages.push('no valid primer pairs found'); return; }
      if (parsed.rows.length > MAX_PAIRS) { messages.push(`too many pairs (${parsed.rows.length}); the maximum is ${MAX_PAIRS}`); return; }
      const pairResults = new Array(parsed.rows.length);
      await runChunked(parsed.rows.length, i => { pairResults[i] = DD.computePair(parsed.rows[i].fwd, parsed.rows[i].rev); }, 'primer pairs');
      const { entries } = DD.pairEntries(parsed.rows, type);   // reuses cached results
      const nc = DD.countNotComputed(pairResults);
      if (nc) messages.push(`${nc} of ${pairResults.length} pairs not computed: ${DD.NOT_COMPUTED}`);
      state.run = { mode: 'pairs', type, rows: parsed.rows, hasOutcome: parsed.hasOutcome, pairResults, entries, stem: state.stems.pairs };
    } else {
      const parsed = DD.parsePrimerList(rows);
      messages.push(...parsed.errors);
      messages.push(...parsed.notes);
      if (parsed.primers.length < 1 || (type !== 'homodimers' && parsed.primers.length < 2)) { messages.push('need at least 2 primers'); return; }
      if (parsed.primers.length > MAX_PRIMERS) { messages.push(`too many primers (${parsed.primers.length}); the maximum is ${MAX_PRIMERS}`); return; }
      const jobs = DD.multiplexJobs(parsed.primers, type);
      const entries = new Array(jobs.length);
      await runChunked(jobs.length, i => { entries[i] = DD.runJob(jobs[i], i + 1); }, 'combinations');
      const nc = DD.countNotComputed(entries.map(e => e.result));
      if (nc) messages.push(`${nc} of ${entries.length} combinations not computed: ${DD.NOT_COMPUTED}`);
      state.run ={ mode: 'multiplex', type, primers: parsed.primers, notes: parsed.notes, entries, stem: state.stems.multiplex };
    }
    render();
    const r = state.run;
    const size = text === EXAMPLES[r.mode] ? 'example' : (r.mode === 'pairs' ? r.rows : r.primers).length;
    track(`run/${r.mode}/${r.type}`, `size/${r.mode}/${size}`);
  } finally {
    state.busy = false;
    $('run-btn').disabled = false;
    setMessages(messages);
  }
}

// ── Rendering ─────────────────────────────────────────────────────────────

function sortedEntries() {
  return DD.sortEntries(state.run.entries, $('sort-mode').value);
}

function alignmentHtml(e, gap) {
  const r = e.result;
  const d = r.which === 'A'
    ? [e.idA, e.idB, e.idB, e.idA] : [e.idB, e.idA, e.idA, e.idB];
  const block = (aln, top, bot) => {
    if (!aln) return `<div class="align-title">3' ${esc(top)} vs ${esc(bot)}</div><span class="muted">no significant structure</span>`;
    const lines = [aln.top_line, aln.match_line, aln.bot_line].map(l => withGap(l, gap));
    return `<div class="align-title">3' ${esc(top)} vs ${esc(bot)}: ${aln.wc_pairs} WC pairs, ${aln.consec_3p} at 3' end</div>` +
      `<div class="align-block">${esc(lines.join('\n'))}</div>`;
  };
  if (r.min_end_hetero === null) return `<span class="muted">not computed: ${esc(DD.notComputedReason(e.seqA, e.seqB))}</span>`;
  let html = block(r.aln, d[0], d[1]);
  if (e.type === 'heterodimer') html += `<details class="alt"><summary>other direction</summary>${block(r.aln_alt, d[2], d[3])}</details>`;
  return html;
}

function render() {
  const runData = state.run;
  if (!runData) return;
  const threshold = getThreshold();
  const gap = $('gap-char').value;
  const all = sortedEntries();
  const flaggedOnly = $('filter-mode').value === 'flagged';
  const shown = flaggedOnly ? all.filter(e => DD.isFlagged(e.result.min_end_hetero, threshold)) : all;

  const counts = { high: 0, medium: 0, low: 0, none: 0 };
  let nFlagged = 0;
  all.forEach(e => {
    counts[DD.riskTier(e.result.min_end_hetero)]++;
    if (DD.isFlagged(e.result.min_end_hetero, threshold)) nFlagged++;
  });
  const what = runData.mode === 'pairs' ? `${runData.rows.length} pairs` : `${runData.primers.length} primers`;
  $('summary-bar').innerHTML =
    `<span>${what}, <strong>${all.length}</strong> combinations</span>` +
    `<span>flagged (&lt; ${esc(fmt2(threshold))} kcal/mol): <strong style="color:#c0392b">${nFlagged}</strong></span>` +
    `<span>${riskBadge(-9)} <strong>${counts.high}</strong></span>` +
    `<span>${riskBadge(-1)} <strong>${counts.medium}</strong></span>` +
    `<span>${riskBadge(0)} <strong>${counts.low}</strong></span>` +
    (counts.none ? `<span>${riskBadge(null)} <strong>${counts.none}</strong></span>` : '') +
    `<span>shown: <strong>${shown.length}</strong></span>`;

  const showType = runData.type === 'both';
  const showOutcome = runData.mode === 'pairs' && runData.hasOutcome;
  $('results-thead').innerHTML = '<tr><th>#</th>' + (showType ? '<th>type</th>' : '') +
    '<th>ID (left)</th><th>sequence (left)</th><th>ID (right)</th><th>sequence (right)</th>' +
    '<th>&Delta;G left&rarr;right</th><th>&Delta;G right&rarr;left</th><th>min &Delta;G<br>(kcal/mol)</th>' +
    '<th>risk</th><th>flagged</th>' + (showOutcome ? '<th>known</th>' : '') + '</tr>';
  const nCols = 10 + (showType ? 1 : 0) + (showOutcome ? 1 : 0);

  const rowsHtml = shown.slice(0, DISPLAY_LIMIT).map(e => {
    const r = e.result;
    const tier = DD.riskTier(r.min_end_hetero);
    const flagged = DD.isFlagged(r.min_end_hetero, threshold);
    const homo = e.type === 'homodimer';
    return '<tr class="data-row">' +
      `<td>${e.n}</td>` + (showType ? `<td>${e.type}</td>` : '') +
      `<td>${esc(e.idA)}</td><td class="seq-cell" title="${esc(e.seqA)}">${esc(e.seqA)}</td>` +
      `<td>${esc(e.idB)}</td><td class="seq-cell" title="${esc(e.seqB)}">${esc(e.seqB)}</td>` +
      `<td class="num">${fmt2(r['end_left-right_dg'])}</td>` +
      `<td class="num">${homo ? '' : fmt2(r['end_right-left_dg'])}</td>` +
      `<td class="num dg-${tier}">${fmt2(r.min_end_hetero)}</td>` +
      `<td>${riskBadge(r.min_end_hetero)}</td><td>${tier === 'none' ? '' : flagged ? 'YES' : 'NO'}</td>` +
      (showOutcome ? `<td>${e.outcome ? esc(e.outcome) : ''}</td>` : '') + '</tr>' +
      `<tr class="aln-row"><td></td><td colspan="${nCols - 1}">${alignmentHtml(e, gap)}</td></tr>`;
  });
  $('results-tbody').innerHTML = rowsHtml.join('');
  $('more-rows').classList.toggle('hidden', shown.length <= DISPLAY_LIMIT);
  $('more-rows').textContent = `showing the first ${DISPLAY_LIMIT} of ${shown.length} rows; the downloads contain all rows.`;

  const multi = runData.mode === 'multiplex';
  $('result-tabs').classList.toggle('hidden', !multi);
  $('dl-summary').classList.toggle('hidden', !multi);
  if (multi) renderSummary(threshold);
  else switchResultTab('pairs');

  const wasHidden = $('results-section').classList.contains('hidden');
  $('results-section').classList.remove('hidden');
  if (wasHidden) $('results-section').scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function primerSummary(threshold) {
  const flagged = state.run.entries.filter(e => e.type === 'heterodimer' && DD.isFlagged(e.result.min_end_hetero, threshold));
  const map = new Map();
  flagged.forEach(e => {
    [[e.idA, e.seqA], [e.idB, e.seqB]].forEach(([id, seq]) => {
      if (!map.has(id)) map.set(id, { id, seq, count: 0, dgSum: 0, dgMin: Infinity });
      const m = map.get(id);
      m.count++; m.dgSum += e.result.min_end_hetero; m.dgMin = Math.min(m.dgMin, e.result.min_end_hetero);
    });
  });
  const rows = [...map.values()].sort((a, b) => b.count - a.count || a.dgSum / a.count - b.dgSum / b.count);
  return { rows, total: flagged.length };
}

function renderSummary(threshold) {
  const { rows, total } = primerSummary(threshold);
  if (!rows.length) {
    $('summary-thead').innerHTML = '';
    $('summary-tbody').innerHTML = '<tr><td style="text-align:center;color:#888;padding:16px">no flagged heterodimers at this threshold</td></tr>';
    return;
  }
  $('summary-thead').innerHTML = '<tr><th>rank</th><th>ID</th><th>sequence</th><th>flagged dimers</th>' +
    '<th>share of flagged (%)</th><th>mean &Delta;G (kcal/mol)</th><th>lowest &Delta;G (kcal/mol)</th></tr>';
  $('summary-tbody').innerHTML = rows.map((r, i) =>
    `<tr><td>${i + 1}</td><td>${esc(r.id)}</td><td class="seq-cell" title="${esc(r.seq)}">${esc(r.seq)}</td>` +
    `<td class="num"><strong>${r.count}</strong></td><td class="num">${(r.count / total * 100).toFixed(1)}</td>` +
    `<td class="num">${fmt2(r.dgSum / r.count)}</td><td class="num">${fmt2(r.dgMin)}</td></tr>`).join('');
}

function switchResultTab(tab) {
  document.querySelectorAll('.result-tab').forEach(b => b.classList.toggle('active', b.dataset.rtab === tab));
  $('rpane-pairs').classList.toggle('hidden', tab !== 'pairs');
  $('rpane-summary').classList.toggle('hidden', tab !== 'summary');
}

// ── Downloads (same layout as the Python script) ──────────────────────────

function resultTable() {
  const t = getThreshold();
  const r = state.run;
  return r.mode === 'pairs'
    ? DD.pairsTable(r.rows, r.pairResults, t, r.hasOutcome)
    : DD.multiplexTable(sortedEntries(), t);
}

function downloadTable(fmt) {
  const table = resultTable();
  if (fmt === 'xlsx') {
    if (typeof XLSX === 'undefined') { setMessages(['Excel export needs SheetJS, which could not be loaded; use .tsv or .csv']); return; }
    const wb = XLSX.utils.book_new();
    XLSX.utils.book_append_sheet(wb, XLSX.utils.aoa_to_sheet(DD.tableToAoa(table)), 'DimerDetective');
    XLSX.writeFile(wb, stem() + '_results.xlsx');
  } else {
    download(DD.tableToText(table, fmt === 'csv' ? ',' : '\t'), stem() + '_results.' + fmt);
  }
}

function downloadAlignments() {
  const r = state.run, t = getThreshold();
  const text = r.mode === 'pairs'
    ? DD.pairsReport(r.rows, r.pairResults, t, r.type)
    : DD.multiplexReport(r.primers.length, sortedEntries(), t, r.type, $('sort-mode').value, r.notes);
  download(reportWithGap(text, $('gap-char').value), stem() + '_alignments.txt');
}

function downloadSummary() {
  const { rows, total } = primerSummary(getThreshold());
  const lines = ['rank\tid\tsequence\tflagged_dimers\tshare_of_flagged_pct\tmean_dg\tlowest_dg\tversion'];
  rows.forEach((r, i) => lines.push([i + 1, r.id, r.seq, r.count, (r.count / total * 100).toFixed(1),
    (r.dgSum / r.count).toFixed(4), r.dgMin.toFixed(4), DD.VERSION_LABEL].join('\t')));
  download(lines.join('\n') + '\n', stem() + '_primer_summary.tsv');
}

// ── Contact form ──────────────────────────────────────────────────────────
// Sent to the developer by email through Web3Forms. The sender's address is
// the reply-to. A hidden field catches bots: if it is filled in, nothing is
// sent but the form reports success.

function showInfoSection(id) {
  document.querySelector('.page-tab[data-page="info"]').click();
  $(id).scrollIntoView({ behavior: 'smooth', block: 'start' });
}

function contactStatus(text, kind) {
  const s = $('contact-status');
  s.textContent = text;
  s.className = 'status' + (kind ? ' ' + kind : '');
}

async function sendContact(e) {
  e.preventDefault();
  const email = $('contact-email').value.trim();
  const subject = $('contact-subject').value.trim();
  const message = $('contact-message').value.trim();
  if (!/^[^\s@]+@[^\s@]+\.[^\s@]+$/.test(email)) { contactStatus('please enter a valid email address', 'err'); return; }
  if (!subject || !message) { contactStatus('please fill in a subject and a message', 'err'); return; }
  if ($('contact-website').value) { $('contact-form').reset(); contactStatus('thank you, your message was sent', 'ok'); return; }
  if (!WEB3FORMS_KEY) { contactStatus('the contact form is not available yet', 'err'); return; }

  $('contact-send').disabled = true;
  contactStatus('sending…');
  try {
    const r = await fetch('https://api.web3forms.com/submit', {
      method: 'POST',
      headers: { 'Content-Type': 'application/json', Accept: 'application/json' },
      body: JSON.stringify({
        access_key: WEB3FORMS_KEY,
        subject: 'DimerDetective: ' + subject,
        from_name: 'DimerDetective contact form',
        email,
        replyto: email,
        message,
        botcheck: false,
      }),
    });
    const data = await r.json().catch(() => ({}));
    if (!r.ok || data.success === false) throw new Error(data.message || 'HTTP ' + r.status);
    $('contact-form').reset();
    contactStatus('thank you, your message was sent', 'ok');
    track('contact/sent');
  } catch (err) {
    contactStatus('sending failed (' + err.message + '); please try again later', 'err');
  } finally {
    $('contact-send').disabled = false;
  }
}

// ── Wiring ────────────────────────────────────────────────────────────────

function wire() {
  document.querySelectorAll('.page-tab').forEach(b => b.addEventListener('click', () => {
    document.querySelectorAll('.page-tab').forEach(x => x.classList.toggle('active', x === b));
    $('page-calc').classList.toggle('hidden', b.dataset.page !== 'calc');
    $('page-info').classList.toggle('hidden', b.dataset.page !== 'info');
    if (b.dataset.page === 'info') track('page/info');
  }));
  document.querySelectorAll('a[href="#contact"]').forEach(a => a.addEventListener('click', e => {
    e.preventDefault();
    showInfoSection('contact');
  }));
  $('contact-form').addEventListener('submit', sendContact);
  if (location.hash === '#contact') setTimeout(() => showInfoSection('contact'), 0);
  document.querySelectorAll('[data-toggle]').forEach(h => h.addEventListener('click', () => {
    const body = $(h.dataset.toggle + '-body');
    const hide = !body.classList.contains('hidden');
    body.classList.toggle('hidden', hide);
    h.querySelector('span').innerHTML = hide ? '&#9654;' : '&#9660;';
  }));
  document.querySelectorAll('.tab[data-mode]').forEach(b => b.addEventListener('click', () => switchMode(b.dataset.mode)));
  document.querySelectorAll('.result-tab').forEach(b => b.addEventListener('click', () => switchResultTab(b.dataset.rtab)));

  $('threshold-preset').addEventListener('change', () => {
    const custom = $('threshold-preset').value === 'custom';
    $('threshold-custom').classList.toggle('hidden', !custom);
    render();
  });
  $('threshold-custom').addEventListener('input', render);
  ['sort-mode', 'filter-mode', 'gap-char'].forEach(id => $(id).addEventListener('change', render));
  const updateHomodimerNote = () => $('homodimer-note').classList.toggle('hidden', $('type-mode').value === 'heterodimers');
  $('type-mode').addEventListener('change', () => {
    updateHomodimerNote();
    if (state.run) setMessages(['interactions changed: run again to update the results']);
  });
  updateHomodimerNote();

  $('run-btn').addEventListener('click', run);
  $('example-btn').addEventListener('click', () => {
    $('input-text').value = EXAMPLES[state.mode];
    state.stems[state.mode] = null;
    setMessages([]);
  });
  $('clear-btn').addEventListener('click', () => {
    $('input-text').value = '';
    state.stems[state.mode] = null;
    setMessages([]);
  });
  $('input-text').addEventListener('keydown', e => {
    if (e.key === 'Enter' && (e.ctrlKey || e.metaKey)) run();
  });

  $('file-input').addEventListener('change', e => { loadFile(e.target.files[0]); e.target.value = ''; });
  const area = $('upload-area');
  area.addEventListener('dragover', e => { e.preventDefault(); area.classList.add('drag'); });
  area.addEventListener('dragleave', () => area.classList.remove('drag'));
  area.addEventListener('drop', e => { e.preventDefault(); area.classList.remove('drag'); loadFile(e.dataTransfer.files[0]); });

  document.querySelectorAll('[data-menu]').forEach(b => b.addEventListener('click', () => $(b.dataset.menu).classList.toggle('open')));
  document.addEventListener('click', e => {
    if (!e.target.closest('.dl-wrap')) document.querySelectorAll('.dl-menu').forEach(m => m.classList.remove('open'));
  });
  document.querySelectorAll('[data-table]').forEach(b => b.addEventListener('click', () => {
    downloadTable(b.dataset.table);
    b.closest('.dl-menu').classList.remove('open');
    track('download/results-' + b.dataset.table);
  }));
  $('dl-alignments').addEventListener('click', () => { downloadAlignments(); track('download/alignments'); });
  $('dl-summary').addEventListener('click', () => { downloadSummary(); track('download/summary'); });

  document.querySelectorAll('.max-pairs').forEach(s => { s.textContent = MAX_PAIRS; });
  document.querySelectorAll('.max-primers').forEach(s => { s.textContent = MAX_PRIMERS; });
  document.querySelectorAll('#version, .version').forEach(s => { s.textContent = DD.VERSION; });
  document.querySelectorAll('.version-label').forEach(s => { s.textContent = DD.VERSION_LABEL; });
  switchMode('pairs');
}

wire();
setupCounting();
DD.init(createThal).then(() => {
  $('engine-status').textContent = 'primer3 engine ready';
  $('run-btn').disabled = false;
}).catch(err => {
  $('engine-status').textContent = 'the primer3 engine could not be loaded: ' + err.message;
});
