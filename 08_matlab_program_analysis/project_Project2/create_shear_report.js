'use strict';
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, HeadingLevel, BorderStyle, WidthType, ShadingType,
  VerticalAlign, ImageRun, PageBreak, Header, Footer, LevelFormat,
  PageNumber, NumberFormat, Tab, TabStopType, TabStopPosition
} = require('docx');
const fs   = require('fs');
const path = require('path');

// ── Paths ────────────────────────────────────────────────────────────────────
const OUT_DOCX = path.join(__dirname, 'ShearDesign_Report.docx');
const IMG_PATH = path.join(__dirname, 'output', 'ShearDesign.png');

// ── Font helpers ─────────────────────────────────────────────────────────────
const F  = 'Times New Roman';
const FH = 'Arial';
const SZ = 24;           // 12 pt body
const SZsm = 22;         // 11 pt for tables
const SZh1 = 28;
const SZh2 = 26;
const SZh3 = 24;

// Convenience: normal run
const r  = (t, o = {}) => new TextRun({ text: t, font: F, size: SZ, color: '000000', ...o });
const rb = (t, o = {}) => r(t, { bold: true, ...o });
const ri = (t, o = {}) => r(t, { italics: true, ...o });
const rsub = (t, o = {}) => r(t, { subScript: true, ...o });
const rsup = (t, o = {}) => r(t, { superScript: true, ...o });
const rbsub = (t, o = {}) => r(t, { bold: true, subScript: true, ...o });

// Table-sized run
const rt  = (t, o = {}) => new TextRun({ text: t, font: F, size: SZsm, color: '000000', ...o });
const rtb = (t, o = {}) => rt(t, { bold: true, ...o });
const rtsub = (t, o = {}) => rt(t, { subScript: true, ...o });

// Greek / math symbols
const LAMBDA = '\u03BB';
const PHI    = '\u03C6';
const ETA    = '\u03B7';
const SQRT   = '\u221A';
const TIMES  = '\u00D7';
const GEQ    = '\u2265';
const LEQ    = '\u2264';
const PRIME  = '\u2032';

// ── Paragraph helpers ────────────────────────────────────────────────────────
function bodyPara(runs, opts = {}) {
  const children = typeof runs === 'string' ? [r(runs)] : runs;
  return new Paragraph({
    children,
    alignment: AlignmentType.BOTH,
    spacing: { before: 60, after: 60, line: 276 },
    ...opts
  });
}

function indentPara(runs, opts = {}) {
  return bodyPara(runs, { indent: { left: 720 }, ...opts });
}

function centeredPara(runs, opts = {}) {
  const children = typeof runs === 'string' ? [r(runs)] : runs;
  return new Paragraph({
    children,
    alignment: AlignmentType.CENTER,
    spacing: { before: 60, after: 60, line: 276 },
    ...opts
  });
}

function blankLine() {
  return new Paragraph({ text: '', spacing: { before: 0, after: 0 } });
}

// Heading helpers — use Arial Bold, black only
function h1(text) {
  return new Paragraph({
    children: [new TextRun({ text, font: FH, size: SZh1, bold: true, color: '000000' })],
    spacing: { before: 360, after: 120 },
    alignment: AlignmentType.LEFT
  });
}
function h2(text) {
  return new Paragraph({
    children: [new TextRun({ text, font: FH, size: SZh2, bold: true, color: '000000' })],
    spacing: { before: 280, after: 100 },
    alignment: AlignmentType.LEFT
  });
}
function h3(text) {
  return new Paragraph({
    children: [new TextRun({ text, font: FH, size: SZh3, bold: true, color: '000000' })],
    spacing: { before: 200, after: 80 },
    alignment: AlignmentType.LEFT
  });
}

// ── Table helpers ────────────────────────────────────────────────────────────
const BD  = { style: BorderStyle.SINGLE, size: 1, color: '000000' };
const BDS = { top: BD, bottom: BD, left: BD, right: BD };
const CELL_MARGINS = { top: 80, bottom: 80, left: 120, right: 120 };
const SHADE_HDR = { fill: 'D9D9D9', type: ShadingType.CLEAR, color: '000000' };
const SHADE_ALT = { fill: 'EEEEEE', type: ShadingType.CLEAR, color: '000000' };

function cell(runs, width, opts = {}) {
  const children = typeof runs === 'string' ? [rt(runs)] : runs;
  const align = opts.center ? AlignmentType.CENTER : (opts.right ? AlignmentType.RIGHT : AlignmentType.LEFT);
  return new TableCell({
    borders: BDS,
    width: { size: width, type: WidthType.DXA },
    margins: CELL_MARGINS,
    shading: opts.shade || undefined,
    verticalAlign: VerticalAlign.CENTER,
    children: [new Paragraph({
      children,
      spacing: { before: 20, after: 20 },
      alignment: align
    })]
  });
}

function hdrCell(runs, width) {
  return cell(runs, width, { shade: SHADE_HDR, center: true });
}

function makeTable(columnWidths, headerTexts, dataRows, opts = {}) {
  const totalWidth = columnWidths.reduce((a, b) => a + b, 0);
  const headerRow = new TableRow({
    children: headerTexts.map((t, i) => hdrCell(
      typeof t === 'string' ? [rtb(t)] : t,
      columnWidths[i]
    )),
    tableHeader: true
  });
  const rows = dataRows.map((row, ri) => new TableRow({
    children: row.map((c, ci) => {
      const shade = opts.altRows && ri % 2 === 1 ? SHADE_ALT : undefined;
      const isFirst = ci === 0;
      return cell(
        typeof c === 'string' ? [rt(c)] : c,
        columnWidths[ci],
        { shade, center: !isFirst && !opts.leftAll, right: false }
      );
    })
  }));
  return new Table({
    columnWidths,
    width: { size: totalWidth, type: WidthType.DXA },
    rows: [headerRow, ...rows]
  });
}

// ── Build subscripted text helpers ───────────────────────────────────────────
// Returns array of TextRun for things like "f'c", "Vci", "bw", "dp" etc.
function fprime_c(sz) {
  const s = sz || SZ;
  return [
    new TextRun({ text: 'f', font: F, size: s, italics: true, color: '000000' }),
    new TextRun({ text: '\u2032', font: F, size: s, color: '000000' }),
    new TextRun({ text: 'c', font: F, size: s, subScript: true, color: '000000' })
  ];
}
function fprime_ci(sz) {
  const s = sz || SZ;
  return [
    new TextRun({ text: 'f', font: F, size: s, italics: true, color: '000000' }),
    new TextRun({ text: '\u2032', font: F, size: s, color: '000000' }),
    new TextRun({ text: 'ci', font: F, size: s, subScript: true, color: '000000' })
  ];
}
function sub(base, subscript, sz, bold) {
  const s = sz || SZ;
  return [
    new TextRun({ text: base, font: F, size: s, color: '000000', bold: bold || false }),
    new TextRun({ text: subscript, font: F, size: s, subScript: true, color: '000000', bold: bold || false })
  ];
}

// Table-sized sub
function tsub(base, subscript, bold) {
  return sub(base, subscript, SZsm, bold);
}

// ── COVER PAGE ───────────────────────────────────────────────────────────────
function coverPage() {
  const ctr = AlignmentType.CENTER;
  const sp = { before: 0, after: 0, line: 276 };
  const ln = (text, bold, size) => new Paragraph({
    alignment: ctr,
    spacing: sp,
    children: [new TextRun({ text, font: F, size: size || 28, bold: bold || false, color: '000000' })]
  });
  return [
    blankLine(), blankLine(), blankLine(), blankLine(), blankLine(),
    blankLine(), blankLine(), blankLine(), blankLine(), blankLine(),
    new Paragraph({
      alignment: ctr, spacing: { before: 0, after: 120 },
      children: [new TextRun({ text: 'Project 1 \u2014 Precast Prestressed Beam', font: F, size: 36, bold: true, color: '000000' })]
    }),
    new Paragraph({
      alignment: ctr, spacing: { before: 0, after: 200 },
      children: [new TextRun({ text: 'Assignment 4: Shear Design', font: F, size: 32, bold: true, color: '000000' })]
    }),
    blankLine(), blankLine(),
    ln('CEE 530 Prestressed Concrete', false, 28),
    ln('Spring 2026', false, 28),
    blankLine(), blankLine(),
    ln('Due: 3/23/2026', false, 26),
    blankLine(), blankLine(), blankLine(),
    ln('Group 3', true, 28),
    ln('Hailey Crosson', false, 26),
    ln('Chidchanok Pleesudjai (Fen)', false, 26),
    blankLine(),
    new Paragraph({ children: [new PageBreak()] })
  ];
}

// ── SECTION 1: Given Data ────────────────────────────────────────────────────
function section1() {
  const content = [];
  content.push(h1('1. Given Data'));

  // 1.1 Section Properties
  content.push(h2('1.1 Section Properties'));
  const W1 = 4600, W2 = 4760;
  const secCols = [W1, W2];
  const secData = [
    [[rt('h')], [rt('28.00 in.')]],
    [[...tsub('A', 'c')], [rt('487.000 in\u00B2')]],
    [[...tsub('I', 'c')], [rt('34,638.8 in\u2074')]],
    [[...tsub('y', 'c'), rt(' (from bottom)')], [rt('20.362 in.')]],
    [[...tsub('y', 'b')], [rt('20.362 in.')]],
    [[...tsub('y', 't')], [rt('7.638 in.')]],
    [[...tsub('b', 'w'), rt(' (sum of stem widths)')], [rt('7.50 in.')]],
    [[...tsub('S', 't'), rt(' = '), ...tsub('I', 'c'), rt('/'), ...tsub('y', 't')], [rt('4,535.1 in\u00B3')]],
    [[...tsub('S', 'b'), rt(' = '), ...tsub('I', 'c'), rt('/'), ...tsub('y', 'b')], [rt('1,701.1 in\u00B3')]],
  ];
  content.push(makeTable(secCols, ['Property', 'Value'], secData, { altRows: true }));

  // 1.2 Material Properties
  content.push(h2('1.2 Material Properties'));
  const matCols = [4600, 4760];
  const matData = [
    [[...fprime_c(SZsm)], [rt('6.0 ksi')]],
    [[...fprime_ci(SZsm)], [rt('4.8 ksi')]],
    [[...tsub('E', 'c')], [rt('4,700 ksi')]],
    [[...tsub('f', 'pu')], [rt('270 ksi')]],
    [[...tsub('f', 'py')], [rt('243 ksi')]],
    [[...tsub('E', 'ps')], [rt('28,500 ksi')]],
    [[...tsub('f', 'y'), rt(' (stirrups)')], [rt('60 ksi')]],
    [[rt(LAMBDA)], [rt('1.0 (normal weight)')]],
    [[rt(PHI), new TextRun({ text: 'v', font: F, size: SZsm, subScript: true, color: '000000' })], [rt('0.75')]],
  ];
  content.push(makeTable(matCols, ['Property', 'Value'], matData, { altRows: true }));

  // 1.3 Prestress Data
  content.push(h2('1.3 Prestress Data'));
  content.push(bodyPara([
    r('Four 0.5-in. diameter Grade 270 strands, each with '),
    ...sub('A', 'ps'), r(' = 0.153 in\u00B2. Initial force per strand '),
    ...sub('P', 'i'), r(' = 25.0 kip; total initial prestress '),
    ...sub('F', 'i'), r(' = 100.0 kip. With 15% losses ('),
    r(ETA), r(' = 0.85), effective prestress '),
    ...sub('P', 'e'), r(' = 85.0 kip.')
  ]));
  content.push(bodyPara([
    r('Tendons 1 and 3 are straight at y = 6.00 in. Tendons 2 and 4 are trapezoidal harped: y = 20.36 in. at supports, sloping to y = 6.00 in. at x = 240 in., flat to x = 528 in., then sloping back to y = 20.36 in.')
  ]));
  content.push(bodyPara([
    r('The vertical component of the effective prestress force, '),
    ...sub('V', 'p'), r(' = '), ...sub('P', 'e'),
    r(' sin(\u03B1), where \u03B1 = arctan(\u0394y/\u0394x) is the tendon inclination angle. For harped tendons 2 and 4 in the sloped region (0 \u2264 x \u2264 240 in.):')
  ]));
  content.push(indentPara([
    r('tan(\u03B1) = (20.36 \u2212 6.00) / 240 = 14.36 / 240 = 0.05983')
  ]));
  content.push(indentPara([
    ...sub('V', 'p'),
    r(' = 2 tendons ' + TIMES + ' 21.25 kip ' + TIMES + ' sin(\u03B1) = 2 ' + TIMES + ' 21.25 ' + TIMES + ' 0.05972 = 2.538 kip  (upward)')
  ]));
  content.push(bodyPara([
    r('In the flat region (240 \u2264 x \u2264 528 in.), the tendon slope is zero and '),
    ...sub('V', 'p'),
    r(' = 0. Thus '),
    ...sub('V', 'p'),
    r(' = 2.538 kip at Sections 1 and 2, and '),
    ...sub('V', 'p'),
    r(' = 0 at Sections 3 and 4.')
  ]));

  // 1.4 Applied Loads
  content.push(h2('1.4 Applied Loads'));
  const ldCols = [4300, 2530, 2530];
  const ldData = [
    [[...tsub('w', 'sw'), rt(' (self-weight)')], [rt('0.04227')], [rt('0.5073')]],
    [[...tsub('w', 'SDL'), rt(' (2-in. topping)')], [rt('0.02083')], [rt('0.2500')]],
    [[...tsub('w', 'LL'), rt(' (live load)')], [rt('0.03500')], [rt('0.4200')]],
    [[...tsub('w', 'DL'), rt(' = '), ...tsub('w', 'sw'), rt(' + '), ...tsub('w', 'SDL')], [rt('0.06311')], [rt('0.7573')]],
    [[...tsub('w', 'u'), rt(' = 1.2'), rt('\u00B7'), ...tsub('w', 'DL'), rt(' + 1.6'), rt('\u00B7'), ...tsub('w', 'LL')], [rt('0.13173')], [rt('1.5808')]],
  ];
  content.push(makeTable(ldCols, ['Load', 'kip/in', 'kip/ft'], ldData, { altRows: true }));

  return content;
}

// ── SECTION 2: Shear Design Method ──────────────────────────────────────────
function section2() {
  const content = [];
  content.push(h1('2. Shear Design Method (CEE 530 / ACI 318-19)'));

  // Introductory explanation
  content.push(bodyPara([
    r('ACI 318-19 requires that the factored shear demand '),
    ...sub('V', 'u'),
    r(' at each section shall not exceed the design shear capacity '),
    r(PHI), ...sub('V', 'n'),
    r(', where '), ...sub('V', 'n'), r(' = '), ...sub('V', 'c'),
    r(' + '), ...sub('V', 's'),
    r('. The concrete contribution '), ...sub('V', 'c'),
    r(' is taken as the lesser of two failure modes: '),
    rb('flexural-shear cracking'), r(' ('), ...sub('V', 'ci'),
    r(') and '), rb('web-shear cracking'), r(' ('), ...sub('V', 'cw'), r(').'),
  ]));

  content.push(bodyPara([
    rb('Flexural-shear cracking ('), ...sub('V', 'ci'), rb(')'),
    r(' initiates as a flexural crack that subsequently turns into a diagonal shear crack. It governs in regions of high moment (near midspan). The computation requires the cracking moment '),
    ...sub('M', 'cr'),
    r(', the unfactored dead-load shear '),
    ...sub('V', 'd'),
    r(', and the factored external-load shear/moment pair '),
    ...sub('V', 'i'), r('/'), ...sub('M', 'max'),
    r('. A lower bound '), ...sub('V', 'ci,min'),
    r(' is enforced to prevent unconservative results when '),
    ...sub('M', 'cr'), r(' is small or zero.'),
  ]));

  content.push(bodyPara([
    rb('Web-shear cracking ('), ...sub('V', 'cw'), rb(')'),
    r(' initiates as a diagonal crack at the neutral axis of the beam before any flexural cracking occurs. It governs near the supports where shear is high and moment is low. The calculation depends on the compressive stress at the centroid '),
    ...sub('f', 'pc'),
    r(' and the vertical component of the prestress force '),
    ...sub('V', 'p'),
    r(' (which acts upward and helps resist applied shear).'),
  ]));

  content.push(bodyPara([
    r('If '), ...sub('V', 'u'), r(' > '), r(PHI), ...sub('V', 'c'),
    r(', transverse reinforcement (stirrups) shall be provided to carry the excess shear '),
    ...sub('V', 's'), r(' = '), ...sub('V', 'u'), r('/'), r(PHI),
    r(' \u2212 '), ...sub('V', 'c'),
    r('. Even when '), ...sub('V', 'u'), r(' \u2264 '), r(PHI), ...sub('V', 'c'),
    r(', ACI 318-19 requires a minimum amount of shear reinforcement ('),
    ...sub('A', 'v'), r('/s)'), ...sub('', 'min'),
    r(' per Table 9.6.3.3.'),
  ]));

  content.push(bodyPara([
    r('The CEE 530 course convention uses '),
    ...sub('V', 'ci,min'),
    r(` = 5${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(`${TIMES}`), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(' (instead of the ACI 318-19 code value of 1.7).'),
  ]));

  // 2.1 Key Equations
  content.push(h2('2.1 Key Equations'));

  // Eq. (1) with explanation
  content.push(bodyPara([
    rb('Flexural-shear strength'), r(' (ACI 318-19 Eq. 22.5.8.3.1a):'),
  ]));
  content.push(indentPara([
    rb('Eq. (1): '),
    ...sub('V', 'ci'), r(` = 0.6${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(TIMES), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(' + '), ...sub('V', 'd'), r(' + ('), ...sub('V', 'i'), r('/'), ...sub('M', 'max'),
    r(`)${TIMES}`), ...sub('M', 'cr'),
    r(` ${GEQ} `), ...sub('V', 'ci,min'),
    r(` = 5${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(TIMES), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
  ]));
  content.push(bodyPara([
    r('where '), ...sub('V', 'd'), r(' = shear from unfactored dead load (SW + SDL); '),
    ...sub('V', 'i'), r(' and '), ...sub('M', 'max'),
    r(' = factored shear and moment from superimposed loads only (1.2'), r(TIMES),
    ...sub('w', 'SDL'), r(' + 1.6'), r(TIMES), ...sub('w', 'LL'), r(').'),
  ]));

  // Eq. (2) with explanation
  content.push(bodyPara([
    rb('Cracking moment'), r(' (ACI 318-19 Eq. 22.5.8.3.2a):'),
  ]));
  content.push(indentPara([
    rb('Eq. (2): '),
    ...sub('M', 'cr'), r(' = ('), ...sub('I', 'c'), r('/'), ...sub('y', 'b'),
    r(`)(6${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(' + '), ...sub('f', 'ce'), r(' \u2212 '), ...sub('f', 'd'), r(')'),
  ]));
  content.push(bodyPara([
    r('where '), ...sub('f', 'ce'), r(' = compressive stress at the bottom fiber due to effective prestress alone ('),
    ...sub('f', 'ce'), r(' = '), ...sub('P', 'e'), r('/'), ...sub('A', 'c'),
    r(' + '), ...sub('P', 'e'), r(TIMES), r('e'), r(TIMES), ...sub('y', 'b'),
    r('/'), ...sub('I', 'c'), r('); '),
    ...sub('f', 'd'), r(' = tensile stress at the bottom fiber due to unfactored dead-load moment ('),
    ...sub('f', 'd'), r(' = '), ...sub('M', 'd'), r(TIMES), ...sub('y', 'b'),
    r('/'), ...sub('I', 'c'), r('). '),
    r('When '), ...sub('f', 'd'), r(' > 6'), r(LAMBDA), r(SQRT), ...fprime_c(),
    r(' + '), ...sub('f', 'ce'),
    r(', the section is already cracked under dead load and '),
    ...sub('M', 'cr'), r(' = 0; in that case '), ...sub('V', 'ci,min'), r(' governs.'),
  ]));

  // Eq. (3) with explanation
  content.push(bodyPara([
    rb('Web-shear strength'), r(' (ACI 318-19 Eq. 22.5.8.3.2):'),
  ]));
  content.push(indentPara([
    rb('Eq. (3): '),
    ...sub('V', 'cw'), r(` = (3.5${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(' + 0.3'), r(TIMES), ...sub('f', 'pc'), r(')'),
    r(TIMES), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(' + '), ...sub('V', 'p'),
  ]));
  content.push(bodyPara([
    r('where '), ...sub('f', 'pc'), r(' = compressive stress at the centroid of the section due to effective prestress ('),
    ...sub('f', 'pc'), r(' = '), ...sub('P', 'e'), r('/'), ...sub('A', 'c'),
    r('); '), ...sub('V', 'p'), r(' = vertical component of the effective prestress force from inclined (harped) tendons.'),
  ]));

  // Eq. (4)
  content.push(bodyPara([
    rb('Concrete shear capacity:'),
  ]));
  content.push(indentPara([
    rb('Eq. (4): '),
    ...sub('V', 'c'), r(' = min('), ...sub('V', 'ci'), r(', '), ...sub('V', 'cw'), r(')'),
  ]));

  // Eq. (5) with explanation
  content.push(bodyPara([
    rb('Required steel shear contribution'), r(' (when '), ...sub('V', 'u'),
    r(' > '), r(PHI), ...sub('V', 'c'), r('):'),
  ]));
  content.push(indentPara([
    rb('Eq. (5): '),
    ...sub('V', 's,req'), r(' = max('), ...sub('V', 'u'), r('/'), r(PHI),
    rsub('v'), r(' \u2212 '), ...sub('V', 'c'), r(', 0)'),
  ]));

  // Eq. (6) with explanation
  content.push(bodyPara([
    rb('Stirrup spacing'), r(' for vertical stirrups (\u03B1'), rsub('s'),
    r(' = 90\u00B0):'),
  ]));
  content.push(indentPara([
    rb('Eq. (6): '),
    ...sub('A', 'v'), r('/s = '), ...sub('V', 's'), r('/('),
    ...sub('f', 'y'), r(TIMES), ...sub('d', 'p'), r(')'),
  ]));
  content.push(bodyPara([
    r('where '), ...sub('A', 'v'), r(' = total cross-sectional area of the stirrup legs; s = stirrup spacing; '),
    ...sub('f', 'y'), r(' = yield strength of stirrup steel; '),
    ...sub('d', 'p'), r(' = effective depth to the prestress centroid.'),
  ]));

  // 2.2 Effective Depth
  content.push(h2('2.2 Effective Depth and Critical Section'));
  content.push(bodyPara([
    r('Per ACI 318-19 Section 22.5.3.2, the effective depth for shear shall be taken as:'),
  ]));
  content.push(indentPara([
    ...sub('d', 'p'),
    r(' = max(h \u2212 '), ...sub('y', 'ps'),
    r(', 0.80'), r(TIMES), r('h) = max(28.00 \u2212 '), ...sub('y', 'ps'),
    r(', 22.40) in.'),
  ]));
  content.push(bodyPara([
    r('Since the tendon centroid '), ...sub('y', 'ps'),
    r(' is relatively high in the section, 0.80'), r(TIMES),
    r('h = 22.40 in. governs at all sections. The critical section for shear is located at '),
    ...sub('d', 'p'), r(' = 22.40 in. (1.87 ft) from the face of the support (ACI 318-19 R22.5.3.2).'),
  ]));

  return content;
}

// ── SECTION 3: Detailed Calculations ─────────────────────────────────────────
function sectionCalc(num, x_in, x_ft, data) {
  // data: { yps, e, fce_term1, fce_term2, fce, fd, Md, fr_fce_fd, mcr,
  //          wext, Vi, Mmax, Vd, term1, term2_vd, term3, vci_sum, vci_min, vci,
  //          fpc, vcw_bracket, Vp, vcw,
  //          Vc, phiVc, Vu, stirrup_needed, Vs, avs_demand, avs_min, s_calc, smax, s_final }
  const content = [];
  content.push(h2(`3.${num} Section ${num}: x = ${x_in} in. (${x_ft} ft)`));

  // Prestress at section
  content.push(h3('Prestress at section'));
  content.push(bodyPara([
    ...sub('y', 'ps'), r(` = ${data.yps} in.,   e = ${data.e} in.`)
  ]));

  // Cracking moment
  content.push(h3('Cracking Moment'));
  content.push(bodyPara([
    ...sub('f', 'ce'), r(' = '), ...sub('P', 'e'), r('/'), ...sub('A', 'c'),
    r(' + '), ...sub('P', 'e'), r(TIMES), r('e'), r(TIMES), ...sub('y', 'b'),
    r('/'), ...sub('I', 'c')
  ]));
  content.push(indentPara([
    r(`= 85.0/487.0 + 85.0${TIMES}${data.e}${TIMES}20.362/34,638.8`),
  ]));
  content.push(indentPara([
    r(`= ${data.fce_term1} + ${data.fce_term2} = ${data.fce} ksi`)
  ]));

  content.push(bodyPara([
    ...sub('f', 'd'), r(' = '), ...sub('M', 'd'), r(TIMES), ...sub('y', 'b'),
    r('/'), ...sub('I', 'c'),
  ]));
  content.push(indentPara([
    r(`where `), ...sub('M', 'd'), r(` = `), ...sub('w', 'DL'),
    r(`${TIMES}x${TIMES}(L\u2212x)/2 = ${data.Md_calc}`)
  ]));
  content.push(indentPara([
    ...sub('f', 'd'), r(` = ${data.fd_calc} = ${data.fd} ksi`)
  ]));

  // fr + fce - fd check
  content.push(bodyPara([
    ...sub('f', 'r'), r(' + '), ...sub('f', 'ce'), r(' \u2212 '), ...sub('f', 'd'),
    r(` = 0.4648 + ${data.fce} \u2212 ${data.fd} = ${data.fr_fce_fd} ksi`),
    r(data.fr_fce_fd_positive
      ? ' > 0'
      : ' < 0')
  ]));

  if (data.fr_fce_fd_positive) {
    content.push(bodyPara([
      ...sub('M', 'cr'), r(' = ('), ...sub('I', 'c'), r('/'), ...sub('y', 'b'),
      r(`)${TIMES}${data.fr_fce_fd} = 1701.1${TIMES}${data.fr_fce_fd} = ${data.mcr} kip-in`)
    ]));
  } else {
    content.push(bodyPara([
      ...sub('M', 'cr'), r(' = 0 (section cracked under dead load alone)')
    ]));
  }

  // Vci
  content.push(h3('Flexural-Shear Strength (Vci)'));

  // Show the main equation FIRST
  content.push(bodyPara([
    r('From Eq. (1):'),
  ]));
  content.push(indentPara([
    ...sub('V', 'ci'), r(` = 0.6${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(TIMES), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(' + '), ...sub('V', 'd'), r(' + ('), ...sub('V', 'i'), r('/'), ...sub('M', 'max'),
    r(`)${TIMES}`), ...sub('M', 'cr'),
  ]));

  // Now compute each component with clear labels
  content.push(bodyPara([
    rb('Factored external-load shear and moment'), r(' (superimposed loads only, excluding self-weight):'),
  ]));
  if (data.wext_show) {
    content.push(indentPara([
      ...sub('w', 'ext'), r(` = 1.2${TIMES}`), ...sub('w', 'SDL'),
      r(` + 1.6${TIMES}`), ...sub('w', 'LL'),
      r(` = 1.2${TIMES}0.02083 + 1.6${TIMES}0.03500 = ${data.wext} kip/in`),
    ]));
  }
  content.push(indentPara([
    ...sub('V', 'i'), r(` = `), ...sub('w', 'ext'),
    r(`${TIMES}(L/2 \u2212 x) = ${data.wext}${TIMES}(384 \u2212 ${x_in}) = ${data.Vi} kip`),
  ]));
  if (data.Mmax_show) {
    content.push(indentPara([
      ...sub('M', 'max'), r(` = `), ...sub('w', 'ext'),
      r(`${TIMES}x${TIMES}(L\u2212x)/2 = ${data.Mmax_calc}`),
    ]));
  }

  content.push(bodyPara([
    rb('Unfactored dead-load shear'), r(' (self-weight + SDL, '), ri('no load factor applied'), r('):'),
  ]));
  content.push(indentPara([
    ...sub('V', 'd'), r(` = `), ...sub('w', 'DL'),
    r(`${TIMES}(L/2 \u2212 x) = 0.06311${TIMES}${data.halfL_minus_x} = ${data.Vd} kip`),
  ]));

  content.push(bodyPara([
    rb('Substitution into Eq. (1):'),
  ]));
  content.push(indentPara([
    ...sub('V', 'ci'),
    r(` = 0.6${TIMES}${LAMBDA}${TIMES}${SQRT}f\u2032`), rsub('c'),
    r(`${TIMES}`), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(` + `), ...sub('V', 'd'),
    r(` + (`), ...sub('V', 'i'), r('/'), ...sub('M', 'max'), r(`)${TIMES}`), ...sub('M', 'cr'),
  ]));
  content.push(indentPara([
    r(`= 0.6${TIMES}1.0${TIMES}0.07746${TIMES}7.50${TIMES}22.40 + ${data.Vd} + ${data.term3}`),
  ]));
  content.push(indentPara([
    r(`= ${data.term1} + ${data.Vd} + ${data.term3} = `),
    rb(`${data.vci_sum} kip`),
  ]));

  content.push(bodyPara([
    rb('Lower bound check:'),
  ]));
  content.push(indentPara([
    ...sub('V', 'ci,min'),
    r(` = 5${TIMES}${LAMBDA}${TIMES}${SQRT}f\u2032`), rsub('c'),
    r(`${TIMES}`), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(` = 5${TIMES}1.0${TIMES}0.07746${TIMES}7.50${TIMES}22.40 = ${data.vci_min} kip`),
  ]));
  content.push(indentPara([
    ...sub('V', 'ci'),
    r(` = max(${data.vci_sum}, ${data.vci_min}) = `),
    rb(`${data.vci} kip`),
  ]));

  // Vcw
  content.push(h3('Vcw calculation'));
  content.push(bodyPara([
    ...sub('f', 'pc'), r(` = `), ...sub('P', 'e'), r('/'), ...sub('A', 'c'),
    r(` = 85.0/487.0 = ${data.fpc} ksi`)
  ]));
  content.push(indentPara([
    r(`(3.5${TIMES}0.07746 + 0.3${TIMES}${data.fpc})${TIMES}7.50${TIMES}22.40 = ${data.vcw_bracket} kip`)
  ]));
  content.push(indentPara([
    ...sub('V', 'p'), r(` = ${data.Vp} kip`)
  ]));
  content.push(bodyPara([
    ...sub('V', 'cw'), r(` = ${data.vcw_bracket} + ${data.Vp} = ${data.vcw} kip`)
  ]));

  // Vc and stirrup check
  content.push(h3('Concrete Shear Capacity and Stirrup Check'));

  // Explain why Vc = min
  content.push(bodyPara([
    r('Per Eq. (4), the concrete shear capacity is taken as the '),
    ri('lesser'), r(' of the two failure modes, because the section will crack by whichever mechanism requires the smaller shear force:'),
  ]));

  // Determine which governs
  const vcw_governs = parseFloat(data.Vc) <= parseFloat(data.vcw) + 0.01;
  const govLabel = vcw_governs ? 'Vcw governs (web-shear cracking controls)' : 'Vci governs (flexural-shear cracking controls)';

  content.push(indentPara([
    ...sub('V', 'c'), r(` = min(`), ...sub('V', 'ci'), r(`, `), ...sub('V', 'cw'),
    r(`) = min(${data.vci}, ${data.vcw}) = `), rb(`${data.Vc} kip`),
  ]));
  content.push(indentPara([
    r(`\u2192 `), rb(govLabel),
  ]));

  content.push(bodyPara([
    r('Design shear capacity:'),
  ]));
  content.push(indentPara([
    r(PHI), ...sub('V', 'c'), r(` = 0.75${TIMES}${data.Vc} = ${data.phiVc} kip`),
  ]));

  content.push(bodyPara([
    r('Factored shear demand at this section:'),
  ]));
  content.push(indentPara([
    ...sub('V', 'u'), r(` = `), ...sub('w', 'u'),
    r(`${TIMES}(L/2 \u2212 x) = 0.13173${TIMES}${data.halfL_minus_x} = ${data.Vu} kip`),
  ]));

  if (data.stirrup_needed) {
    content.push(bodyPara([
      ...sub('V', 'u'), r(` = ${data.Vu} kip > `), r(PHI), ...sub('V', 'c'),
      r(` = ${data.phiVc} kip \u2192 `), rb('Shear reinforcement required beyond minimum'),
    ]));

    content.push(bodyPara([
      r('Required steel contribution from Eq. (5):'),
    ]));
    content.push(indentPara([
      ...sub('V', 's'), r(` = `), ...sub('V', 'u'), r(`/`), r(PHI), rsub('v'),
      r(` \u2212 `), ...sub('V', 'c'),
      r(` = ${data.Vu}/0.75 \u2212 ${data.Vc} = `), rb(`${data.Vs} kip`),
    ]));

    content.push(bodyPara([
      r('Stirrup sizing: use '), rb('#3 closed U-stirrups'),
      r(` (2 legs ${TIMES} 0.11 in\u00B2/leg = `), ...sub('A', 'v'), r(' = 0.22 in\u00B2):'),
    ]));
    content.push(bodyPara([
      r('Required '), ...sub('A', 'v'), r('/s from demand (Eq. 6):'),
    ]));
    content.push(indentPara([
      ...sub('A', 'v'), r(`/s`), rsub('demand'), r(` = `), ...sub('V', 's'),
      r(`/(`), ...sub('f', 'y'), r(TIMES), ...sub('d', 'p'),
      r(`) = ${data.Vs}/(60${TIMES}22.40) = ${data.avs_demand} in\u00B2/in`),
    ]));
    content.push(bodyPara([
      r('Minimum '), ...sub('A', 'v'), r('/s (from Section 4, Criterion 1 governs):'),
    ]));
    content.push(indentPara([
      ...sub('A', 'v'), r(`/s`), rsub('min'), r(` = 0.00726 in\u00B2/in > ${data.avs_demand} in\u00B2/in \u2192 `),
      rb('minimum governs'),
    ]));
    content.push(bodyPara([
      r('Spacing:'),
    ]));
    content.push(indentPara([
      r(`s`), rsub('req'), r(` = `), ...sub('A', 'v'), r(`/(A`), rsub('v'), r(`/s)`),
      r(` = 0.22/0.00726 = ${data.s_calc} in.`),
    ]));
    content.push(bodyPara([
      r('Maximum spacing per ACI 318-19 Section 9.7.6.2.2:'),
    ]));
    content.push(indentPara([
      ...sub('s', 'max'), r(` = min(3h/4, 24 in.) = min(3${TIMES}28/4, 24) = min(21.0, 24) = `),
      rb('21.0 in.'),
    ]));
    content.push(indentPara([
      r('s = min('), rsub('req'), r(`, `), ...sub('s', 'max'),
      r(`) = min(${data.s_calc}, ${data.smax}) = `),
      rb(`${data.s_final} in.`), r(' \u2190 '), rb('USE'),
    ]));
  } else {
    content.push(bodyPara([
      ...sub('V', 'u'), r(` = ${data.Vu} kip ${LEQ} `), r(PHI), ...sub('V', 'c'),
      r(` = ${data.phiVc} kip \u2192 `), rb('Minimum stirrups sufficient'),
    ]));
    content.push(bodyPara([
      r('Use '), rb('#3 closed U-stirrups'),
      r(` (2 legs ${TIMES} 0.11 in\u00B2/leg = `), ...sub('A', 'v'),
      r(' = 0.22 in\u00B2). Spacing governed by minimum '), ...sub('A', 'v'),
      r('/s and ACI 318-19 Section 9.7.6.2.2:'),
    ]));
    content.push(indentPara([
      ...sub('s', 'max'), r(` = min(3h/4, 24 in.) = min(21.0, 24) = 21.0 in.`),
    ]));
    content.push(indentPara([
      r('s = min(0.22/0.00726, 21.0) = min(30.30, 21.0) = '),
      rb(`${data.s_final} in.`), r(' \u2190 '), rb('USE'),
    ]));
  }

  return content;
}

function section3() {
  const content = [];
  content.push(h1('3. Detailed Calculations at Design Sections'));

  // Section 1: x = 18 in.
  content.push(...sectionCalc(1, 18, '1.5', {
    yps: '12.642', e: '7.721',
    fce_term1: '0.1745', fce_term2: '0.3858', fce: '0.5603',
    Md_calc: `0.06311${TIMES}18${TIMES}(768\u221218)/2 = 426.0 kip-in`,
    fd_calc: `426.0${TIMES}20.362/34,638.8`,
    fd: '0.2504',
    fr_fce_fd: '0.7747', fr_fce_fd_positive: true,
    mcr: '1,317.8',
    wext_show: true, wext: '0.08100',
    Vi: '29.646',
    Mmax_show: true,
    Mmax_calc: `0.08100${TIMES}18${TIMES}750/2 = 546.75 kip-in`,
    halfL_minus_x: '366',
    Vd: '23.097',
    term1: '7.808', term3: '71.454',
    vci_sum: '102.36', vci_min: '65.07', vci: '102.36',
    fpc: '0.1745',
    vcw_bracket: '54.343', Vp: '2.538', vcw: '56.881',
    Vc: '56.88', phiVc: '42.66',
    Vu: '48.21',
    stirrup_needed: true,
    Vs: '7.40', avs_demand: '0.00551', avs_min: '0.00726',
    s_calc: '30.30', smax: '21.0', s_final: '21'
  }));

  // Section 2: x = 120 in.
  content.push(...sectionCalc(2, 120, '10.0', {
    yps: '9.590', e: '10.772',
    fce_term1: '0.1745', fce_term2: '0.5383', fce: '0.7128',
    Md_calc: `0.06311${TIMES}120${TIMES}(768\u2212120)/2 = 2,453.5 kip-in`,
    fd_calc: `2,453.5${TIMES}20.362/34,638.8`,
    fd: '1.4423',
    fr_fce_fd: '\u22120.2648', fr_fce_fd_positive: false,
    mcr: '0',
    wext_show: false, wext: '0.08100',
    Vi: '21.384',
    Mmax_show: false,
    Mmax_calc: '',
    halfL_minus_x: '264',
    Vd: '16.660',
    term1: '7.808', term3: '0',
    vci_sum: '24.47', vci_min: '65.07', vci: '65.07',
    fpc: '0.1745',
    vcw_bracket: '54.343', Vp: '2.538', vcw: '56.881',
    Vc: '56.88', phiVc: '42.66',
    Vu: '34.78',
    stirrup_needed: false,
    Vs: '0', avs_demand: '0', avs_min: '0.00726',
    s_calc: '', smax: '21.0', s_final: '21'
  }));

  // Section 3: x = 240 in.
  content.push(...sectionCalc(3, 240, '20.0', {
    yps: '6.000', e: '14.362',
    fce_term1: '0.1745', fce_term2: '0.7177', fce: '0.8922',
    Md_calc: `0.06311${TIMES}240${TIMES}(768\u2212240)/2 = 3,997.4 kip-in`,
    fd_calc: `3,997.4${TIMES}20.362/34,638.8`,
    fd: '2.3505',
    fr_fce_fd: '\u22120.9936', fr_fce_fd_positive: false,
    mcr: '0',
    wext_show: false, wext: '0.08100',
    Vi: '11.664',
    Mmax_show: false,
    Mmax_calc: '',
    halfL_minus_x: '144',
    Vd: '9.088',
    term1: '7.808', term3: '0',
    vci_sum: '16.90', vci_min: '65.07', vci: '65.07',
    fpc: '0.1745',
    vcw_bracket: '54.343', Vp: '0.000', vcw: '54.343',
    Vc: '54.34', phiVc: '40.76',
    Vu: '18.97',
    stirrup_needed: false,
    Vs: '0', avs_demand: '0', avs_min: '0.00726',
    s_calc: '', smax: '21.0', s_final: '21'
  }));

  // Section 4: x = 384 in.
  content.push(...sectionCalc(4, 384, '32.0', {
    yps: '6.000', e: '14.362',
    fce_term1: '0.1745', fce_term2: '0.7177', fce: '0.8922',
    Md_calc: `0.06311${TIMES}384${TIMES}(768\u2212384)/2 = 4,651.1 kip-in`,
    fd_calc: `4,651.1${TIMES}20.362/34,638.8`,
    fd: '2.7351',
    fr_fce_fd: '\u22121.3781', fr_fce_fd_positive: false,
    mcr: '0',
    wext_show: false, wext: '0.08100',
    Vi: '0.000',
    Mmax_show: false,
    Mmax_calc: '',
    halfL_minus_x: '0',
    Vd: '0.000',
    term1: '7.808', term3: '0',
    vci_sum: '7.81', vci_min: '65.07', vci: '65.07',
    fpc: '0.1745',
    vcw_bracket: '54.343', Vp: '0.000', vcw: '54.343',
    Vc: '54.34', phiVc: '40.76',
    Vu: '0.00',
    stirrup_needed: false,
    Vs: '0', avs_demand: '0', avs_min: '0.00726',
    s_calc: '', smax: '21.0', s_final: '21'
  }));

  return content;
}

// ── SECTION 4: Minimum Shear Reinforcement ──────────────────────────────────
function section4() {
  const content = [];
  content.push(h1('4. Minimum Shear Reinforcement (ACI 318-19 Table 9.6.3.3)'));

  // Criterion 1
  content.push(bodyPara([rb('Criterion 1:')]));
  content.push(indentPara([
    ...sub('A', 'v'), r(`/s = 0.75${TIMES}${SQRT}`),
    ...fprime_c(), r(`${TIMES}`), ...sub('b', 'w'), r('/'), ...sub('f', 'y')
  ]));
  content.push(indentPara([
    r(`= 0.75${TIMES}0.07746${TIMES}7.50/60 = 0.00726 in\u00B2/in`),
    r('  '), rb('\u2190 GOVERNS')
  ]));

  // Criterion 2
  content.push(bodyPara([rb('Criterion 2:')]));
  content.push(indentPara([
    ...sub('A', 'v'), r(`/s = (50/1000)${TIMES}`), ...sub('b', 'w'), r('/'), ...sub('f', 'y')
  ]));
  content.push(indentPara([
    r(`= 0.050${TIMES}7.50/60 = 0.00625 in\u00B2/in`)
  ]));

  // Criterion 3
  content.push(bodyPara([rb('Criterion 3:')]));
  content.push(indentPara([
    ...sub('A', 'v'), r(`/s = (`), ...sub('A', 'ps'), r(TIMES), ...sub('f', 'pu'),
    r(`)/(80${TIMES}`), ...sub('f', 'y'), r(TIMES), ...sub('d', 'p'),
    r(`)${TIMES}${SQRT}(`), ...sub('d', 'p'), r('/'), ...sub('b', 'w'), r(')')
  ]));
  content.push(indentPara([
    r(`= (0.612${TIMES}270)/(80${TIMES}60${TIMES}22.4)${TIMES}${SQRT}(22.4/7.5) = 0.00266 in\u00B2/in`)
  ]));

  content.push(blankLine());
  content.push(bodyPara([
    r('Selected bar: #3 closed U-stirrups, '),
    ...sub('A', 'v'), r(` = 2${TIMES}0.11 = 0.22 in\u00B2`)
  ]));

  content.push(blankLine());
  content.push(bodyPara([rb('Spacing limits:')]));
  content.push(indentPara([
    r('Basic: '), ...sub('s', 'max'),
    r(` = min(3h/4, 24) = min(21.0, 24) = 21.0 in.`)
  ]));
  content.push(indentPara([
    r('Tight (when '), ...sub('V', 's'), r(` > 4${TIMES}${SQRT}`),
    ...fprime_c(), r(TIMES), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r('): '), ...sub('s', 'max'), r(' = min(3h/8, 12) = 10.5 in. \u2014 not triggered')
  ]));

  return content;
}

// ── SECTION 5: Summary ──────────────────────────────────────────────────────
function section5() {
  const content = [];
  content.push(h1('5. Summary and Design'));

  // Summary table
  const cols = [900, 900, 1100, 1100, 1100, 1000, 1000, 1360, 900];
  const hdr = [
    [rtb('Section')],
    [rtb('x (ft)')],
    [rtb('V'), new TextRun({ text: 'u', font: F, size: SZsm, subScript: true, bold: true, color: '000000' }), rtb(' (kip)')],
    [rtb('V'), new TextRun({ text: 'ci', font: F, size: SZsm, subScript: true, bold: true, color: '000000' }), rtb(' (kip)')],
    [rtb('V'), new TextRun({ text: 'cw', font: F, size: SZsm, subScript: true, bold: true, color: '000000' }), rtb(' (kip)')],
    [rtb('V'), new TextRun({ text: 'c', font: F, size: SZsm, subScript: true, bold: true, color: '000000' }), rtb(' (kip)')],
    [rt(PHI), rtb('V'), new TextRun({ text: 'c', font: F, size: SZsm, subScript: true, bold: true, color: '000000' }), rtb(' (kip)')],
    [rtb('Stirrup')],
    [rtb('s (in.)')]
  ];
  const data = [
    ['1', '1.5', '48.21', '102.36', '56.88', '56.88', '42.66', '#3 U @ ', '21'],
    ['2', '10.0', '34.78', '65.07', '56.88', '56.88', '42.66', '#3 U @ ', '21'],
    ['3', '20.0', '18.97', '65.07', '54.34', '54.34', '40.76', '#3 U @ ', '21'],
    ['4', '32.0', '0.00', '65.07', '54.34', '54.34', '40.76', '#3 U @ ', '21'],
  ];
  content.push(makeTable(cols, hdr, data, { altRows: true }));

  content.push(blankLine());

  // Design recommendation
  content.push(h3('Design Recommendation'));
  content.push(bodyPara([
    ...sub('V', 'cw'),
    r(' governs at all four design sections. Only Section 1 (x = 1.5 ft) requires stirrups beyond the minimum ('),
    ...sub('V', 's'), r(' = 7.40 kip), but minimum '),
    ...sub('A', 'v'), r('/s still governs the spacing. The recommended stirrup design is:'),
  ]));
  content.push(indentPara([
    rb('#3 closed U-stirrups'), r(` (2 legs ${TIMES} 0.11 in\u00B2 = `),
    ...sub('A', 'v'), r(' = 0.22 in\u00B2), '),
    ...sub('f', 'y'), r(' = 60 ksi'),
  ]));
  content.push(indentPara([
    rb('Spacing: s = 21 in. throughout the span'),
    r(' (governed by '), ...sub('s', 'max'), r(' = min(3h/4, 24 in.) per ACI 318-19 Sec. 9.7.6.2.2)'),
  ]));
  content.push(bodyPara([
    r('This uniform spacing satisfies all four design sections.')
  ]));

  // Max Vs check
  content.push(h3('Maximum Vs Check'));
  content.push(bodyPara([
    ...sub('V', 's,max'), r(` = 8${TIMES}${LAMBDA}${TIMES}${SQRT}`),
    ...fprime_c(), r(TIMES), ...sub('b', 'w'), r(TIMES), ...sub('d', 'p'),
    r(` = 8${TIMES}1.0${TIMES}0.07746${TIMES}7.50${TIMES}22.40 = 104.1 kip. At Section 1, `),
    ...sub('V', 's'), r(` = 7.40 kip ${LEQ} 104.1 kip. Section size is adequate.`)
  ]));

  return content;
}

// ── SECTION 6: Figure ────────────────────────────────────────────────────────
function section6() {
  const content = [];
  content.push(h1('6. Shear Design Figure'));

  // Insert image
  let imgData;
  try {
    imgData = fs.readFileSync(IMG_PATH);
  } catch (e) {
    console.warn('WARNING: Could not read ShearDesign.png at', IMG_PATH);
    content.push(bodyPara('[Image not found: ShearDesign.png]'));
    return content;
  }

  content.push(new Paragraph({
    alignment: AlignmentType.CENTER,
    spacing: { before: 120, after: 60 },
    children: [
      new ImageRun({
        data: imgData,
        type: 'png',
        transformation: { width: 624, height: 468 }  // ~6.5 in wide
      })
    ]
  }));

  content.push(centeredPara([
    rb('Fig. 1'), r(' \u2014 Shear demand, concrete capacity, required '),
    ...sub('A', 'v'), r('/s, and stirrup spacing along the beam span.')
  ]));

  return content;
}

// ── ACI Code References ─────────────────────────────────────────────────────
function aciReferences() {
  const content = [];
  content.push(h1('ACI Code References'));
  const refs = [
    ['Table 22.5.8.2', 'Detailed method (Vci/Vcw)'],
    ['Eq. 22.5.8.3.1a', 'Vci flexural-shear'],
    ['Eq. 22.5.8.3.2', 'Vcw web-shear'],
    ['Sec. 22.5.3.2', 'Critical section at dp'],
    ['Table 9.6.3.3', 'Minimum Av/s'],
    ['Sec. 9.7.6.2.2', 'Spacing limits'],
    ['Table 21.2.1', '\u03C6 = 0.75 for shear'],
  ];
  const cols = [3400, 5960];
  content.push(makeTable(cols, ['ACI 318-19 Reference', 'Description'], refs, { altRows: true, leftAll: true }));
  return content;
}

// ── BUILD DOCUMENT ──────────────────────────────────────────────────────────
async function main() {
  const children = [
    ...coverPage(),
    ...section1(),
    ...section2(),
    ...section3(),
    ...section4(),
    ...section5(),
    ...section6(),
    ...aciReferences()
  ];

  const doc = new Document({
    styles: {
      default: {
        document: {
          run: { font: F, size: SZ, color: '000000' },
          paragraph: { spacing: { line: 276 } }
        }
      },
      paragraphStyles: [
        {
          id: 'Heading1',
          name: 'Heading 1',
          basedOn: 'Normal',
          next: 'Normal',
          run: { font: FH, size: SZh1, bold: true, color: '000000' },
          paragraph: { spacing: { before: 360, after: 120 } }
        },
        {
          id: 'Heading2',
          name: 'Heading 2',
          basedOn: 'Normal',
          next: 'Normal',
          run: { font: FH, size: SZh2, bold: true, color: '000000' },
          paragraph: { spacing: { before: 280, after: 100 } }
        },
        {
          id: 'Heading3',
          name: 'Heading 3',
          basedOn: 'Normal',
          next: 'Normal',
          run: { font: FH, size: SZh3, bold: true, color: '000000' },
          paragraph: { spacing: { before: 200, after: 80 } }
        }
      ]
    },
    sections: [{
      properties: {
        page: {
          size: { width: 12240, height: 15840, orientation: 'portrait' },
          margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 }
        },
        pageNumberStart: 1,
        pageNumberFormatType: NumberFormat.DECIMAL
      },
      headers: {
        default: new Header({
          children: [new Paragraph({
            alignment: AlignmentType.LEFT,
            children: [
              new TextRun({ text: 'CEE 530 \u2014 Assignment 4: Shear Design', font: F, size: 20, color: '000000', italics: true })
            ]
          })]
        })
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [
              new TextRun({ text: 'Page ', font: F, size: 20, color: '000000' }),
              new TextRun({ children: [PageNumber.CURRENT], font: F, size: 20, color: '000000' })
            ]
          })]
        })
      },
      children
    }]
  });

  const buffer = await Packer.toBuffer(doc);
  fs.writeFileSync(OUT_DOCX, buffer);
  console.log('Report written to:', OUT_DOCX);
  console.log('Size:', (buffer.length / 1024).toFixed(1), 'KB');
}

main().catch(err => { console.error(err); process.exit(1); });
