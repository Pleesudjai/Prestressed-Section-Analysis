// create_shear_report_v2.js — Shear Design Report (calc-report template)
// CEE 530 Prestressed Concrete — Project 2 — Assignment 4
// Run: node create_shear_report_v2.js

const fs = require("fs");
const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  ImageRun, Header, Footer, AlignmentType, HeadingLevel, BorderStyle,
  WidthType, ShadingType, PageBreak, PageNumber, LevelFormat,
  VerticalAlign,
} = require("docx");

// ── Helpers ─────────────────────────────────────────────────────────
const B  = (t, opts={}) => new TextRun({ text: t, bold: true, font: "Arial", size: 22, ...opts });
const R  = (t, opts={}) => new TextRun({ text: t, font: "Arial", size: 22, ...opts });
const Ri = (t, opts={}) => new TextRun({ text: t, font: "Arial", size: 22, italics: true, ...opts });
const sub = (t, opts={}) => new TextRun({ text: t, font: "Arial", size: 22, subScript: true, ...opts });
const sup = (t, opts={}) => new TextRun({ text: t, font: "Arial", size: 22, superScript: true, ...opts });
const mono = (t, opts={}) => new TextRun({ text: t, font: "Courier New", size: 20, ...opts });

const P = (children, opts={}) => new Paragraph({ children, spacing: { after: 80 }, ...opts });
const Pn = (children, opts={}) => new Paragraph({ children, spacing: { after: 0 }, ...opts });
const Ps = (children, opts={}) => new Paragraph({ children, spacing: { before: 120, after: 80 }, ...opts });
const blank = () => new Paragraph({ children: [], spacing: { after: 40 } });

const heading1 = (text) => new Paragraph({
  children: [new TextRun({ text, bold: true, font: "Arial", size: 28, allCaps: true })],
  spacing: { before: 300, after: 160 },
  border: { bottom: { style: BorderStyle.SINGLE, size: 6, color: "000000" } },
});

const heading2 = (text) => new Paragraph({
  children: [new TextRun({ text, bold: true, font: "Arial", size: 24 })],
  spacing: { before: 200, after: 100 },
});

const heading3 = (text) => new Paragraph({
  children: [new TextRun({ text, bold: true, font: "Arial", size: 22 })],
  spacing: { before: 160, after: 80 },
});

const hrLine = () => new Paragraph({
  children: [],
  border: { bottom: { style: BorderStyle.SINGLE, size: 2, color: "999999" } },
  spacing: { after: 80 },
});

// Table helpers
const border = { style: BorderStyle.SINGLE, size: 1, color: "999999" };
const borders = { top: border, bottom: border, left: border, right: border };
const cellP = (text, opts={}) => new Paragraph({
  children: [new TextRun({ text, font: "Arial", size: 20, ...opts })],
  spacing: { after: 0 },
  alignment: opts.align || AlignmentType.LEFT,
});
const cellPr = (runs, opts={}) => new Paragraph({
  children: runs,
  spacing: { after: 0 },
  alignment: opts.align || AlignmentType.LEFT,
});

const hdrShade = { fill: "E0E0E0", type: ShadingType.CLEAR };
const cellMargins = { top: 40, bottom: 40, left: 80, right: 80 };

function makeCell(content, width, opts={}) {
  const children = typeof content === "string" ? [cellP(content, opts)] : content;
  return new TableCell({
    borders,
    width: { size: width, type: WidthType.DXA },
    margins: cellMargins,
    shading: opts.shading || undefined,
    verticalAlign: opts.vAlign || VerticalAlign.CENTER,
    children,
  });
}

function makeRow(cells) { return new TableRow({ children: cells }); }

// Equation line: monospace with indent
const eq = (text) => new Paragraph({
  children: [new TextRun({ text, font: "Courier New", size: 20 })],
  indent: { left: 480 },
  spacing: { after: 40 },
});

// Result line: bold monospace
const eqResult = (text) => new Paragraph({
  children: [new TextRun({ text, font: "Courier New", size: 20, bold: true })],
  indent: { left: 480 },
  spacing: { after: 80 },
});

// PASS/FAIL line
const checkLine = (text, pass) => new Paragraph({
  children: [
    new TextRun({ text: text + "  \u2192  ", font: "Courier New", size: 20 }),
    new TextRun({ text: pass ? "PASS" : "FAIL", font: "Arial", size: 22, bold: true }),
  ],
  indent: { left: 480 },
  spacing: { after: 80 },
});

// ── Images ──────────────────────────────────────────────────────────
const imgDir = "unpacked_fen/word/media/";
const img1 = fs.readFileSync(imgDir + "image1.png"); // cross-section
const img2 = fs.readFileSync(imgDir + "image2.png"); // beam elevation
const img3 = fs.readFileSync(imgDir + "image3.png"); // shear plots

// Also try MATLAB shear figure if available
let imgShear;
try { imgShear = fs.readFileSync("output/ShearDesign.png"); } catch(e) { imgShear = img3; }

// ── BUILD DOCUMENT ──────────────────────────────────────────────────
const doc = new Document({
  styles: {
    default: {
      document: { run: { font: "Arial", size: 22 } },
    },
  },
  numbering: {
    config: [{
      reference: "bullets",
      levels: [{ level: 0, format: LevelFormat.BULLET, text: "\u2022", alignment: AlignmentType.LEFT,
        style: { paragraph: { indent: { left: 720, hanging: 360 } } } }],
    }],
  },
  sections: [
    // ════════════════════════════════════════════════════════════════
    //  COVER PAGE
    // ════════════════════════════════════════════════════════════════
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },
        },
      },
      children: [
        blank(), blank(), blank(), blank(), blank(),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 200 },
          children: [new TextRun({ text: "CEE 530 PRESTRESSED CONCRETE", font: "Arial", size: 36, bold: true })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 200 },
          children: [new TextRun({ text: "Project 1 \u2014 Assignment 4: Shear Design", font: "Arial", size: 30, bold: true })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 100 },
          children: [new TextRun({ text: "Shear Design of Prestressed Double-T Beam", font: "Arial", size: 26 })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 80 },
          children: [new TextRun({ text: "ACI 318-19 / Naaman Section 5.7 \u2014 Empirical Method (Vci / Vcw)", font: "Arial", size: 22 })],
        }),
        blank(), blank(),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 80 },
          children: [new TextRun({ text: "Spring 2026   |   Due: 3/23/2026", font: "Arial", size: 22 })],
        }),
        blank(),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 60 },
          children: [B("Group 3", { size: 24 })],
        }),
        ...["Hailey Crosson", "Reyaneh Kazemian", "Leohuns Robert Kambu",
            "Prince Aminu", "Chidchanok Pleesudjai (Fen)", "Haribhushan Sirigineedi"
        ].map(name => new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 40 },
          children: [R(name)],
        })),
        blank(), blank(),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          children: [R("Code Edition: CEE 530 (ACI 318-19 shear provisions)")],
        }),
      ],
    },

    // ════════════════════════════════════════════════════════════════
    //  MAIN BODY
    // ════════════════════════════════════════════════════════════════
    {
      properties: {
        page: {
          size: { width: 12240, height: 15840 },
          margin: { top: 1440, right: 1296, bottom: 1440, left: 1296 },
        },
      },
      headers: {
        default: new Header({
          children: [new Paragraph({
            alignment: AlignmentType.RIGHT,
            children: [new TextRun({ text: "CEE 530 \u2014 Project 1: Shear Design", font: "Arial", size: 18, color: "666666" })],
          })],
        }),
      },
      footers: {
        default: new Footer({
          children: [new Paragraph({
            alignment: AlignmentType.CENTER,
            children: [
              new TextRun({ text: "Page ", font: "Arial", size: 18, color: "666666" }),
              new TextRun({ children: [PageNumber.CURRENT], font: "Arial", size: 18, color: "666666" }),
            ],
          })],
        }),
      },
      children: [
        // ── SIGN CONVENTION ──────────────────────────────────────
        new Paragraph({
          spacing: { after: 60 },
          border: { top: { style: BorderStyle.SINGLE, size: 2, color: "000000" },
                    bottom: { style: BorderStyle.SINGLE, size: 2, color: "000000" } },
          children: [
            B("SIGN CONVENTION: "),
            R("Compression = positive (+)     Tension = negative (\u2212)"),
          ],
        }),
        P([
          B("UNIT SYSTEM: "), R("kip, in. (moments in kip-in., stresses in ksi)"),
        ]),
        P([
          B("METHOD: "), R("Detailed method \u2014 Vci / Vcw (ACI 318-19 Table 22.5.8.2)"),
        ]),

        // ══════════════════════════════════════════════════════════
        //  PART 1 \u2014 GIVEN DATA
        // ══════════════════════════════════════════════════════════
        heading1("Part 1 \u2014 Given Data"),

        heading2("1.1  Beam Geometry"),
        P([R("Section type: Double-T (precast, simply supported)")]),
        P([R("Span: L = (20 + 24 + 20) ft \u00D7 12 in./ft = "), B("768 in. (64.0 ft)")]),
        P([R("Total depth: h = "), B("28.00 in.")]),

        // Beam elevation figure
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { before: 120, after: 40 },
          children: [new ImageRun({
            type: "png", data: img2,
            transformation: { width: 360, height: 160 },
            altText: { title: "Beam Elevation", description: "Simply supported beam 20+24+20 ft", name: "beam" },
          })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 120 },
          children: [Ri("Fig. 1 \u2014 Beam elevation and support conditions", { size: 20 })],
        }),

        heading2("1.2  Section Properties"),
        // Cross-section figure
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { before: 80, after: 40 },
          children: [new ImageRun({
            type: "png", data: img1,
            transformation: { width: 420, height: 280 },
            altText: { title: "Cross Section", description: "Double-T cross section with dimensions", name: "xsection" },
          })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 120 },
          children: [Ri("Fig. 2 \u2014 Section dimensions and vertex numbering", { size: 20 })],
        }),

        // Section properties table
        new Table({
          width: { size: 5400, type: WidthType.DXA },
          columnWidths: [2700, 2700],
          rows: [
            makeRow([
              makeCell("Property", 2700, { shading: hdrShade, bold: true }),
              makeCell("Value", 2700, { shading: hdrShade, bold: true }),
            ]),
            ...([
              ["Ac", "487.000 in\u00B2"],
              ["Ic", "34,638.82 in\u2074"],
              ["yc (from bottom)", "20.362 in."],
              ["yb", "20.362 in."],
              ["yt", "7.638 in."],
              ["bw (sum of stem widths)", "7.50 in."],
              ["St = Ic/yt", "4,535.1 in\u00B3"],
              ["Sb = Ic/yb", "1,701.1 in\u00B3"],
            ].map(([prop, val]) => makeRow([
              makeCell(prop, 2700),
              makeCell(val, 2700),
            ]))),
          ],
        }),

        heading2("1.3  Material Properties"),
        new Table({
          width: { size: 5400, type: WidthType.DXA },
          columnWidths: [2700, 2700],
          rows: [
            makeRow([
              makeCell("Property", 2700, { shading: hdrShade, bold: true }),
              makeCell("Value", 2700, { shading: hdrShade, bold: true }),
            ]),
            ...([
              ["f\u2032c", "6.0 ksi"],
              ["f\u2032ci", "4.8 ksi"],
              ["Ec", "4,700 ksi"],
              ["fpu", "270 ksi"],
              ["fpy", "243 ksi"],
              ["Eps", "28,500 ksi"],
              ["fy (stirrups)", "60 ksi"],
              ["\u03BB", "1.0 (normal weight)"],
              ["\u03C6v", "0.75 (ACI Table 21.2.1)"],
            ].map(([prop, val]) => makeRow([
              makeCell(prop, 2700),
              makeCell(val, 2700),
            ]))),
          ],
        }),
        P([R("\u221A(f\u2032c in psi) = \u221A(6000) / 1000 = "), B("0.07746 ksi"), R("  (used throughout shear equations)")]),

        heading2("1.4  Prestress Data"),
        P([R("Four 0.5-in. diameter Grade 270 strands, each Aps = 0.153 in\u00B2.")]),
        eq("Aps,total = 4 \u00D7 0.153 = 0.612 in\u00B2"),
        eq("Pi per strand = 25.0 kip   \u2192   Fi,total = 4 \u00D7 25.0 = 100.0 kip"),
        eq("fpi = Pi/Aps = 25.0 / 0.153 = 163.4 ksi"),
        eq("Loss factor: \u03B7 = 1 \u2212 0.15 = 0.85"),
        eqResult("Pe,total = \u03B7 \u00D7 Fi = 0.85 \u00D7 100.0 = 85.0 kip"),
        blank(),
        P([B("Tendon layout (y from bottom):")]),
        P([R("Tendons 1, 3 (straight): y = 6.00 in. (constant)")]),
        P([R("Tendons 2, 4 (harped): y = 20.36 in. at supports, sloping to y = 6.00 in. at x = 240 in., flat to x = 528 in., then sloping back.")]),
        blank(),
        P([B("Vertical prestress component Vp"), R(" (harped tendons, sloped region 0 \u2264 x \u2264 240 in.):")]),
        eq("tan(\u03B1) = (20.36 \u2212 6.00) / 240 = 14.36 / 240 = 0.05983"),
        eq("Pe per harped tendon = 0.153 \u00D7 163.4 \u00D7 0.85 = 21.25 kip"),
        eq("Vp = 2 \u00D7 21.25 \u00D7 sin(\u03B1) = 2 \u00D7 21.25 \u00D7 0.05972"),
        eqResult("Vp = 2.538 kip (upward)        [0 in flat zone, x \u2265 240 in.]"),

        heading2("1.5  Applied Loads"),
        P([B("Self-weight:")]),
        eq("w_sw = Ac \u00D7 (150/1728) / 1000 = 487.0 \u00D7 0.08681 / 1000"),
        eqResult("w_sw = 0.04227 kip/in. (0.5073 kip/ft)"),
        P([B("Superimposed dead load (2-in. topping, 150 pcf):")]),
        eq("w_SDL = 2 \u00D7 120 \u00D7 (150/1728) / 1000"),
        eqResult("w_SDL = 0.02083 kip/in. (0.2500 kip/ft)"),
        P([B("Live load:")]),
        eq("w_LL = 0.42 kip/ft / 12 = 0.03500 kip/in. (0.4200 kip/ft)"),
        blank(),
        P([B("Total unfactored dead load:")]),
        eq("w_DL = w_sw + w_SDL = 0.04227 + 0.02083 = 0.06311 kip/in."),
        P([B("Factored load (ACI 318-19 Sec. 5.3.1b):")]),
        eq("wu = 1.2 \u00D7 w_DL + 1.6 \u00D7 w_LL"),
        eq("   = 1.2 \u00D7 0.06311 + 1.6 \u00D7 0.03500"),
        eq("   = 0.07573 + 0.05600"),
        eqResult("wu = 0.13173 kip/in. (1.5808 kip/ft)"),

        // ══════════════════════════════════════════════════════════
        //  PART 2 \u2014 SHEAR DESIGN METHOD
        // ══════════════════════════════════════════════════════════
        heading1("Part 2 \u2014 Shear Design Method"),
        P([R("ACI 318-19 requires that the factored shear demand Vu at each section shall not exceed the design shear capacity \u03C6Vn, where Vn = Vc + Vs. The concrete contribution Vc is taken as the lesser of two failure modes: "),
           B("flexural-shear cracking (Vci)"), R(" and "), B("web-shear cracking (Vcw)"), R(".")]),
        blank(),
        P([B("Flexural-shear cracking (Vci)"), R(" initiates as a flexural crack that subsequently turns into a diagonal shear crack. It governs in regions of high moment (near midspan). A lower bound Vci,min is enforced.")]),
        blank(),
        P([B("Web-shear cracking (Vcw)"), R(" initiates as a diagonal crack at the neutral axis before any flexural cracking occurs. It governs near the supports where shear is high and moment is low.")]),
        blank(),

        heading2("2.1  Key Equations"),
        P([B("Eq. (1) \u2014 Flexural-shear strength"), R(" (ACI 318-19 Eq. 22.5.8.3.1a):")]),
        eq("Vci = 0.6\u03BB\u221A(f\u2032c) \u00D7 bw \u00D7 dp  +  Vd  +  (Vi/Mmax) \u00D7 Mcr"),
        eq("\u2265 Vci,min = 1.7\u03BB\u221A(f\u2032c) \u00D7 bw \u00D7 dp      [ACI 318-19 hard floor]"),
        eq("\u2264 Vci,max = 5.0\u03BB\u221A(f\u2032c) \u00D7 bw \u00D7 dp      [CEE 530 convention \u2014 governs]"),
        P([R("where Vd = shear from unfactored dead load (SW + SDL); Vi and Mmax = factored shear and moment from superimposed loads only (1.2\u00D7wSDL + 1.6\u00D7wLL).")]),
        blank(),
        P([B("Eq. (2) \u2014 Cracking moment"), R(" (ACI 318-19 Eq. 22.5.8.3.2a):")]),
        eq("Mcr = (Ic/yb)(6\u03BB\u221A(f\u2032c) + fce \u2212 fd)"),
        P([R("where fce = compressive stress at the bottom fiber due to effective prestress alone; fd = tensile stress at the bottom fiber due to unfactored dead-load moment. When fd > 6\u03BB\u221A(f\u2032c) + fce, Mcr = 0 and Vci,min governs.")]),
        blank(),
        P([B("Eq. (3) \u2014 Web-shear strength"), R(" (ACI 318-19 Eq. 22.5.8.3.2):")]),
        eq("Vcw = (3.5\u03BB\u221A(f\u2032c) + 0.3fpc) \u00D7 bw \u00D7 dp  +  Vp"),
        P([R("where fpc = Pe/Ac; Vp = vertical component of harped tendons.")]),
        blank(),
        P([B("Eq. (4):"), R("  Vc = min(Vci, Vcw)")]),
        P([B("Eq. (5):"), R("  Vs,req = max(Vu/\u03C6v \u2212 Vc, 0)")]),
        P([B("Eq. (6):"), R("  Av/s = Vs / (fy \u00D7 dp)")]),

        heading2("2.2  Bounds on Vci"),
        P([B("Lower bound (ACI 318-19 hard floor):"), R("  Prevents Vci from dropping unrealistically low near midspan where Vi \u2192 0 or Mcr = 0:")]),
        eq("Vci \u2265 Vci,min = 1.7\u03BB\u221A(f\u2032c) \u00D7 bw \u00D7 dp"),
        P([R("Without this floor, Vci would reduce to just Vd + 0.6\u03BB\u221A(f\u2032c)\u00D7bw\u00D7dp in regions where DL tension exceeds prestress + rupture (Mcr = 0) \u2014 an unconservatively small value.")]),
        blank(),
        P([B("Upper bound (CEE 530, Prof. Mobasher):"), R("  Caps Vci near supports where Vi/Mmax \u2192 \u221E would inflate the three-term formula:")]),
        eq("Vci \u2264 Vci,max = 5.0\u03BB\u221A(f\u2032c) \u00D7 bw \u00D7 dp"),
        P([R("The 5.0 coefficient equals the upper cap on the simplified method Vc (ACI 318-19 Sec. 22.5.8.1). Using it as the Vci ceiling in the detailed method ensures that Vci can never exceed the simplified-method maximum \u2014 a course convention.")]),
        blank(),
        P([B("For this beam (bw = 7.50 in., dp = 22.40 in., f\u2032c = 6.0 ksi):")]),
        eq("Vci,min (ACI 318-19) = 1.7 \u00D7 1.0 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 22.12 kip  [floor]"),
        eq("Vci,max (CEE 530)    = 5.0 \u00D7 1.0 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 65.07 kip  [cap]"),

        // ══════════════════════════════════════════════════════════
        //  PART 3 \u2014 DETAILED CALCULATIONS AT SECTION 1
        // ══════════════════════════════════════════════════════════
        heading1("Part 3 \u2014 Detailed Calculations at Section 1: x = 18 in. (1.5 ft)"),

        heading2("3.1  Effective Depth dp"),
        P([R("dp shall not be less than 0.80h per ACI 318-19 Sec. 22.5.3.2.")]),
        eq("Tendon centroid at x = 18 in.:"),
        eq("  T1,T3 (straight): y = 6.00 in."),
        eq("  T2,T4 (harped):   y = 20.36 + (6.00\u221220.36)\u00D7(18/240) = 19.28 in."),
        eq("  y_ps = (6.00 + 19.28 + 6.00 + 19.28) / 4 = 12.642 in."),
        eq("  e = yc \u2212 y_ps = 20.362 \u2212 12.642 = 7.721 in."),
        blank(),
        eq("dp = max(h \u2212 y_ps, 0.80h)"),
        eq("   = max(28.00 \u2212 12.642, 0.80 \u00D7 28.00)"),
        eq("   = max(15.358, 22.40)"),
        eqResult("dp = 22.40 in.   [0.80h governs]"),

        heading2("3.2  Modulus of Rupture fr"),
        P([R("Per ACI 318-19 Table 22.5.8.2, fr for Vci/Mcr uses 6\u03BB\u221A(f\u2032c in psi).")]),
        eq("fr = 6 \u00D7 1.0 \u00D7 \u221A(6000) / 1000 = 6 \u00D7 77.46 / 1000"),
        eqResult("fr = 0.4648 ksi"),

        heading2("3.3  Prestress Stress at Bottom Fiber fce"),
        P([R("fce = compressive stress at bottom fiber due to effective prestress alone (positive = compression).")]),
        eq("Formula:  fce = Pe/Ac + Pe \u00D7 e \u00D7 yb / Ic"),
        blank(),
        eq("Substitution:"),
        eq("  fce = 85.0/487.0 + 85.0 \u00D7 7.721 \u00D7 20.362 / 34,638.82"),
        eq("      = 0.1745 + 13,359.3 / 34,638.82"),
        eq("      = 0.1745 + 0.3858"),
        eqResult("fce = 0.5603 ksi (compression)"),

        heading2("3.4  Dead-Load Stress at Bottom Fiber fd"),
        eq("Md = w_DL \u00D7 x \u00D7 (L\u2212x) / 2"),
        eq("   = 0.06311 \u00D7 18 \u00D7 (768\u221218) / 2"),
        eq("   = 0.06311 \u00D7 18 \u00D7 750 / 2"),
        eqResult("Md = 426.0 kip-in."),
        blank(),
        eq("fd = Md \u00D7 yb / Ic = 426.0 \u00D7 20.362 / 34,638.82"),
        eqResult("fd = 0.2504 ksi"),

        heading2("3.5  Cracking Moment Mcr"),
        eq("fr + fce \u2212 fd = 0.4648 + 0.5603 \u2212 0.2504 = 0.7747 ksi  > 0  \u2714"),
        blank(),
        eq("Mcr = (Ic/yb) \u00D7 (fr + fce \u2212 fd)"),
        eq("    = 1,701.1 \u00D7 0.7747"),
        eqResult("Mcr = 1,317.8 kip-in. (109.8 kip-ft)"),

        heading2("3.6  Flexural-Shear Strength Vci"),
        P([R("Load terms at x = 18 in.:")]),
        eq("w_ext = 1.2 \u00D7 w_SDL + 1.6 \u00D7 w_LL = 1.2\u00D70.02083 + 1.6\u00D70.03500 = 0.08100 kip/in."),
        eq("Vi   = w_ext \u00D7 (L/2 \u2212 x) = 0.08100 \u00D7 (384 \u2212 18) = 0.08100 \u00D7 366 = 29.646 kip"),
        eq("Mmax = w_ext \u00D7 x \u00D7 (L\u2212x)/2 = 0.08100 \u00D7 18 \u00D7 750/2 = 546.75 kip-in."),
        eq("Vd   = w_DL \u00D7 (L/2 \u2212 x) = 0.06311 \u00D7 366 = 23.097 kip"),
        blank(),
        P([R("Three-term Vci from Eq. (1):")]),
        eq("Term 1: 0.6 \u00D7 1.0 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 7.808 kip"),
        eq("Term 2: Vd = 23.097 kip"),
        eq("Term 3: (Vi/Mmax) \u00D7 Mcr = (29.646/546.75) \u00D7 1,317.8 = 71.454 kip"),
        eq("Sum:    7.808 + 23.097 + 71.454 = 102.36 kip"),
        blank(),
        P([R("Bounds check:")]),
        eq("Vci,min (ACI floor) = 1.7 \u00D7 1.0 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 22.12 kip"),
        eq("Vci,max (CEE cap)   = 5.0 \u00D7 1.0 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 65.07 kip"),
        eq("Vci (computed)      = 102.36 kip"),
        checkLine("Floor:  102.36 kip  \u2265  22.12 kip", true),
        checkLine("Cap:    102.36 kip  >  65.07 kip  \u2192  capped to Vci,max", false),
        eqResult("Vci = min(102.36, 65.07) = 65.07 kip"),

        heading2("3.7  Web-Shear Strength Vcw"),
        P([B("Centroidal prestress stress fpc:")]),
        eq("fpc = Pe/Ac = 85.0/487.0 = 0.1745 ksi"),
        blank(),
        P([B("Vertical prestress component Vp"), R(" (x = 18 in., in sloped region):")]),
        eq("Vp = 2.538 kip (upward)"),
        blank(),
        P([B("Vcw from Eq. (3):")]),
        eq("(3.5 \u00D7 0.07746 + 0.3 \u00D7 0.1745) \u00D7 7.50 \u00D7 22.40"),
        eq("= (0.27111 + 0.05235) \u00D7 168.00 = 54.343 kip"),
        eq("Vcw = 54.343 + 2.538"),
        eqResult("Vcw = 56.881 kip"),

        heading2("3.8  Concrete Shear Capacity Vc"),
        eq("Vc = min(Vci, Vcw) = min(65.07, 56.881)"),
        eqResult("Vc = 56.88 kip   \u2192 Vcw GOVERNS (web-shear cracking controls)"),
        eq("\u03C6Vc = 0.75 \u00D7 56.88 = 42.66 kip"),

        heading2("3.9  Shear Demand vs. Capacity"),
        eq("Vu = wu \u00D7 (L/2 \u2212 x) = 0.13173 \u00D7 (384 \u2212 18) = 0.13173 \u00D7 366"),
        eqResult("Vu = 48.21 kip"),
        blank(),
        P([B("Code check (ACI 318-19 Sec. 22.5.1.1):")]),
        checkLine("Vu = 48.21 kip  >  \u03C6Vc = 42.66 kip", false),
        eqResult("STATUS: Shear reinforcement required beyond minimum"),

        heading2("3.10  Required Steel Contribution Vs"),
        eq("Vs,req = max(Vu/\u03C6v \u2212 Vc, 0)"),
        eq("       = max(48.21/0.75 \u2212 56.88, 0)"),
        eq("       = max(64.28 \u2212 56.88, 0)"),
        eqResult("Vs,req = 7.40 kip"),
        blank(),
        P([B("Maximum Vs check (ACI 318-19 Sec. 22.5.1.2):")]),
        eq("Vs,max = 8 \u00D7 \u03BB \u00D7 \u221A(f\u2032c) \u00D7 bw \u00D7 dp = 8 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 104.1 kip"),
        checkLine("Vs,req = 7.40 kip  \u2264  Vs,max = 104.1 kip  (section size adequate)", true),

        heading2("3.11  Stirrup Design"),
        P([R("Selected bar: "), B("#3 closed U-stirrups"), R(" (2 legs \u00D7 0.11 in.\u00B2/leg = Av = 0.22 in.\u00B2)")]),
        blank(),
        P([B("Av/s from demand (Eq. 6):")]),
        eq("Av/s_demand = Vs/(fy \u00D7 dp) = 7.40/(60 \u00D7 22.40) = 7.40/1344"),
        eqResult("Av/s_demand = 0.00551 in.\u00B2/in."),
        blank(),
        P([B("Av/s minimum"), R(" (from Part 4, Criterion 1 governs): "), B("0.00726 in.\u00B2/in.")]),
        eq("Av/s_min = 0.00726 > 0.00551  \u2192  minimum governs"),
        blank(),
        P([B("Spacing:")]),
        eq("s_req = Av / (Av/s) = 0.22 / 0.00726 = 30.30 in."),
        eq("s_max = min(3h/4, 24 in.) = min(21.0, 24) = 21.0 in.  [ACI Sec. 9.7.6.2.2]"),
        eq("s = min(30.30, 21.0)"),
        eqResult("s = 21 in.   \u2190 USE"),
        blank(),
        P([B("Verify capacity:")]),
        eq("Vs,prov = Av \u00D7 fy \u00D7 dp / s = 0.22 \u00D7 60 \u00D7 22.40 / 21 = 14.08 kip"),
        eq("\u03C6(Vc + Vs) = 0.75 \u00D7 (56.88 + 14.08) = 0.75 \u00D7 70.96 = 53.22 kip"),
        checkLine("\u03C6(Vc + Vs) = 53.22 kip  \u2265  Vu = 48.21 kip", true),

        // ══════════════════════════════════════════════════════════
        //  DETAILED CALCULATIONS AT SECTION 2
        // ══════════════════════════════════════════════════════════
        new Paragraph({ children: [new PageBreak()] }),
        heading1("Part 3B \u2014 Detailed Calculations at Section 2: x = 120 in. (10.0 ft)"),

        heading2("Effective Depth dp"),
        eq("Tendon centroid at x = 120 in.:"),
        eq("  T1,T3 (straight): y = 6.00 in."),
        eq("  T2,T4 (harped):   y = 20.36 + (6.00\u221220.36)\u00D7(120/240) = 13.18 in."),
        eq("  y_ps = (6.00 + 13.18 + 6.00 + 13.18) / 4 = 9.590 in."),
        eq("  e = yc \u2212 y_ps = 20.362 \u2212 9.590 = 10.772 in."),
        blank(),
        eq("dp = max(h \u2212 y_ps, 0.80h) = max(28.00 \u2212 9.590, 22.40) = max(18.41, 22.40)"),
        eqResult("dp = 22.40 in.   [0.80h governs]"),

        heading2("Prestress & Dead-Load Stresses"),
        eq("fce = Pe/Ac + Pe \u00D7 e \u00D7 yb / Ic"),
        eq("    = 85.0/487.0 + 85.0 \u00D7 10.772 \u00D7 20.362 / 34,638.82"),
        eq("    = 0.1745 + 0.5383"),
        eqResult("fce = 0.7128 ksi"),
        blank(),
        eq("Md = 0.06311 \u00D7 120 \u00D7 (768\u2212120) / 2 = 0.06311 \u00D7 120 \u00D7 648/2"),
        eqResult("Md = 2,453.0 kip-in."),
        eq("fd = 2,453.0 \u00D7 20.362 / 34,638.82"),
        eqResult("fd = 1.4423 ksi"),

        heading2("Cracking Moment Mcr"),
        eq("fr + fce \u2212 fd = 0.4648 + 0.7128 \u2212 1.4423 = \u22120.2648 ksi  < 0  \u2716"),
        eqResult("Mcr = 0   (DL tension exceeds prestress + rupture)"),

        heading2("Flexural-Shear Strength Vci"),
        eq("Vi   = 0.08100 \u00D7 (384 \u2212 120) = 0.08100 \u00D7 264 = 21.384 kip"),
        eq("Mmax = 0.08100 \u00D7 120 \u00D7 648/2 = 3,149.3 kip-in."),
        eq("Vd   = 0.06311 \u00D7 264 = 16.661 kip"),
        blank(),
        eq("Term 1: 0.6 \u00D7 0.07746 \u00D7 7.50 \u00D7 22.40 = 7.808 kip"),
        eq("Term 2: Vd = 16.661 kip"),
        eq("Term 3: (Vi/Mmax) \u00D7 Mcr = 0  (Mcr = 0)"),
        eq("Sum:    7.808 + 16.661 + 0 = 24.47 kip"),
        blank(),
        P([R("Bounds check:")]),
        eq("Vci,min (ACI floor) = 22.12 kip    Vci,max (CEE cap) = 65.07 kip"),
        checkLine("Floor:  24.47 kip  \u2265  22.12 kip", true),
        checkLine("Cap:    24.47 kip  \u2264  65.07 kip", true),
        eqResult("Vci = 24.47 kip   (no bound governs; Mcr = 0 limits the value)"),

        heading2("Web-Shear Strength Vcw"),
        eq("fpc = 0.1745 ksi     Vp = 2.538 kip (x = 120 in., sloped region)"),
        eq("Vcw = (3.5\u00D70.07746 + 0.3\u00D70.1745) \u00D7 7.50 \u00D7 22.40 + 2.538"),
        eq("    = 54.343 + 2.538"),
        eqResult("Vcw = 56.881 kip"),

        heading2("Vc, Vu, Vs, Stirrup"),
        eq("Vc = min(Vci, Vcw) = min(24.47, 56.881)"),
        eqResult("Vc = 24.47 kip   \u2192 Vci GOVERNS (flexural-shear controls)"),
        eq("\u03C6Vc = 0.75 \u00D7 24.47 = 18.35 kip"),
        blank(),
        eq("Vu = 0.13173 \u00D7 (384 \u2212 120) = 0.13173 \u00D7 264"),
        eqResult("Vu = 34.78 kip"),
        checkLine("Vu = 34.78 kip  >  \u03C6Vc = 18.35 kip", false),
        eqResult("STATUS: Shear reinforcement required"),
        blank(),
        eq("Vs,req = max(34.78/0.75 \u2212 24.47, 0) = max(46.37 \u2212 24.47, 0)"),
        eqResult("Vs,req = 21.90 kip"),
        checkLine("Vs,req = 21.90 kip  \u2264  Vs,max = 104.1 kip", true),
        blank(),
        eq("Av/s_demand = 21.90 / (60 \u00D7 22.40) = 21.90/1344"),
        eqResult("Av/s_demand = 0.01629 in.\u00B2/in.   > Av/s_min = 0.00726  \u2192  demand governs"),
        eq("s_req = 0.22 / 0.01629 = 13.50 in."),
        eq("s_max = 21.0 in.   [Vs = 21.90 < 4\u221Af\u2032c\u00D7bw\u00D7dp = 52.06 \u2192 basic limit]"),
        eqResult("s = min(13.50, 21.0) = 13 in.   \u2190 USE"),
        blank(),
        P([B("Verify capacity:")]),
        eq("Vs,prov = 0.22 \u00D7 60 \u00D7 22.40 / 13 = 22.74 kip"),
        eq("\u03C6(Vc + Vs) = 0.75 \u00D7 (24.47 + 22.74) = 0.75 \u00D7 47.21 = 35.41 kip"),
        checkLine("\u03C6(Vc + Vs) = 35.41 kip  \u2265  Vu = 34.78 kip", true),

        // ══════════════════════════════════════════════════════════
        //  DETAILED CALCULATIONS AT SECTION 3
        // ══════════════════════════════════════════════════════════
        new Paragraph({ children: [new PageBreak()] }),
        heading1("Part 3C \u2014 Detailed Calculations at Section 3: x = 240 in. (20.0 ft)"),

        heading2("Effective Depth dp"),
        eq("Tendon centroid at x = 240 in. (drape point):"),
        eq("  T1,T3 (straight): y = 6.00 in."),
        eq("  T2,T4 (harped):   y = 6.00 in. (at drape point)"),
        eq("  y_ps = (6.00 + 6.00 + 6.00 + 6.00) / 4 = 6.000 in."),
        eq("  e = 20.362 \u2212 6.000 = 14.362 in."),
        blank(),
        eq("dp = max(28.00 \u2212 6.00, 22.40) = max(22.00, 22.40)"),
        eqResult("dp = 22.40 in.   [0.80h governs]"),

        heading2("Prestress & Dead-Load Stresses"),
        eq("fce = 85.0/487.0 + 85.0 \u00D7 14.362 \u00D7 20.362 / 34,638.82"),
        eq("    = 0.1745 + 0.7177"),
        eqResult("fce = 0.8922 ksi"),
        blank(),
        eq("Md = 0.06311 \u00D7 240 \u00D7 528/2"),
        eqResult("Md = 3,997.7 kip-in."),
        eq("fd = 3,997.7 \u00D7 20.362 / 34,638.82"),
        eqResult("fd = 2.3505 ksi"),

        heading2("Cracking Moment Mcr"),
        eq("fr + fce \u2212 fd = 0.4648 + 0.8922 \u2212 2.3505 = \u22120.9936 ksi  < 0  \u2716"),
        eqResult("Mcr = 0   (DL tension exceeds prestress + rupture)"),

        heading2("Flexural-Shear Strength Vci"),
        eq("Vi   = 0.08100 \u00D7 (384 \u2212 240) = 0.08100 \u00D7 144 = 11.664 kip"),
        eq("Mmax = 0.08100 \u00D7 240 \u00D7 528/2 = 5,132.2 kip-in."),
        eq("Vd   = 0.06311 \u00D7 144 = 9.088 kip"),
        blank(),
        eq("Term 1: 7.808 kip     Term 2: 9.088 kip     Term 3: 0 (Mcr = 0)"),
        eq("Sum:    7.808 + 9.088 + 0 = 16.90 kip"),
        blank(),
        P([R("Bounds check:")]),
        eq("Vci,min (ACI floor) = 22.12 kip    Vci,max (CEE cap) = 65.07 kip"),
        checkLine("Floor:  16.90 kip  <  22.12 kip  \u2192  raised to Vci,min", false),
        eqResult("Vci = max(16.90, 22.12) = 22.12 kip   (ACI floor governs)"),

        heading2("Web-Shear Strength Vcw"),
        eq("fpc = 0.1745 ksi     Vp = 0 kip (x = 240 in., flat zone)"),
        eq("Vcw = (3.5\u00D70.07746 + 0.3\u00D70.1745) \u00D7 7.50 \u00D7 22.40 + 0"),
        eqResult("Vcw = 54.343 kip"),

        heading2("Vc, Vu, Vs, Stirrup"),
        eq("Vc = min(22.12, 54.343)"),
        eqResult("Vc = 22.12 kip   \u2192 Vci GOVERNS (ACI floor controls)"),
        eq("\u03C6Vc = 0.75 \u00D7 22.12 = 16.59 kip"),
        blank(),
        eq("Vu = 0.13173 \u00D7 (384 \u2212 240) = 0.13173 \u00D7 144"),
        eqResult("Vu = 18.97 kip"),
        checkLine("Vu = 18.97 kip  >  \u03C6Vc = 16.59 kip", false),
        eqResult("STATUS: Shear reinforcement required"),
        blank(),
        eq("Vs,req = max(18.97/0.75 \u2212 22.12, 0) = max(25.29 \u2212 22.12, 0)"),
        eqResult("Vs,req = 3.17 kip"),
        eq("Av/s_demand = 3.17/1344 = 0.00236 in.\u00B2/in."),
        eq("Av/s_min = 0.00726  >  0.00236  \u2192  minimum governs"),
        eq("s_req = 0.22 / 0.00726 = 30.30 in.    s_max = 21.0 in."),
        eqResult("s = min(30.30, 21.0) = 21 in.   \u2190 USE"),
        blank(),
        P([B("Verify capacity:")]),
        eq("Vs,prov = 0.22 \u00D7 60 \u00D7 22.40 / 21 = 14.08 kip"),
        eq("\u03C6(Vc + Vs) = 0.75 \u00D7 (22.12 + 14.08) = 0.75 \u00D7 36.20 = 27.15 kip"),
        checkLine("\u03C6(Vc + Vs) = 27.15 kip  \u2265  Vu = 18.97 kip", true),

        // ══════════════════════════════════════════════════════════
        //  DETAILED CALCULATIONS AT SECTION 4
        // ══════════════════════════════════════════════════════════
        new Paragraph({ children: [new PageBreak()] }),
        heading1("Part 3D \u2014 Detailed Calculations at Section 4: x = 384 in. (32.0 ft, midspan)"),

        heading2("Effective Depth dp"),
        eq("Tendon centroid at x = 384 in. (midspan, flat zone):"),
        eq("  All tendons at y = 6.00 in.   \u2192  y_ps = 6.000 in."),
        eq("  e = 20.362 \u2212 6.000 = 14.362 in."),
        eq("dp = max(22.00, 22.40)"),
        eqResult("dp = 22.40 in.   [0.80h governs]"),

        heading2("Prestress & Dead-Load Stresses"),
        eq("fce = 85.0/487.0 + 85.0 \u00D7 14.362 \u00D7 20.362 / 34,638.82"),
        eqResult("fce = 0.8922 ksi   (same as Section 3)"),
        blank(),
        eq("Md = 0.06311 \u00D7 384 \u00D7 384/2"),
        eqResult("Md = 4,653.2 kip-in."),
        eq("fd = 4,653.2 \u00D7 20.362 / 34,638.82"),
        eqResult("fd = 2.7351 ksi"),

        heading2("Cracking Moment Mcr"),
        eq("fr + fce \u2212 fd = 0.4648 + 0.8922 \u2212 2.7351 = \u22121.3782 ksi  < 0  \u2716"),
        eqResult("Mcr = 0   (DL tension exceeds prestress + rupture)"),

        heading2("Flexural-Shear Strength Vci"),
        eq("Vi   = 0.08100 \u00D7 (384 \u2212 384) = 0 kip   (midspan)"),
        eq("Vd   = 0.06311 \u00D7 (384 \u2212 384) = 0 kip   (midspan)"),
        blank(),
        eq("Term 1: 7.808 kip     Term 2: 0 kip     Term 3: 0 (Mcr = 0)"),
        eq("Sum:    7.808 kip"),
        blank(),
        P([R("Bounds check:")]),
        eq("Vci,min (ACI floor) = 22.12 kip    Vci,max (CEE cap) = 65.07 kip"),
        checkLine("Floor:  7.808 kip  <  22.12 kip  \u2192  raised to Vci,min", false),
        eqResult("Vci = max(7.808, 22.12) = 22.12 kip   (ACI floor governs)"),

        heading2("Web-Shear Strength Vcw"),
        eq("fpc = 0.1745 ksi     Vp = 0 kip (flat zone)"),
        eqResult("Vcw = 54.343 kip   (same as Section 3)"),

        heading2("Vc, Vu, Vs, Stirrup"),
        eq("Vc = min(22.12, 54.343)"),
        eqResult("Vc = 22.12 kip   \u2192 Vci GOVERNS (ACI floor controls)"),
        eq("\u03C6Vc = 0.75 \u00D7 22.12 = 16.59 kip"),
        blank(),
        eqResult("Vu = 0 kip   (midspan of simply supported beam)"),
        checkLine("Vu = 0 kip  \u2264  \u03C6Vc = 16.59 kip", true),
        eqResult("STATUS: Only minimum stirrups required"),
        blank(),
        eq("Vs,req = 0 kip"),
        eq("Av/s_min = 0.00726 governs    s = min(0.22/0.00726, 21.0) = min(30.30, 21.0)"),
        eqResult("s = 21 in.   \u2190 USE"),

        // ══════════════════════════════════════════════════════════
        //  PART 4 \u2014 MINIMUM SHEAR REINFORCEMENT
        // ══════════════════════════════════════════════════════════
        new Paragraph({ children: [new PageBreak()] }),
        heading1("Part 4 \u2014 Minimum Shear Reinforcement (ACI 318-19 Table 9.6.3.3)"),
        P([R("Three criteria shall be evaluated; the largest governs.")]),
        blank(),
        P([B("Criterion 1 (Eq. 9.6.3.3a):")]),
        eq("Av/s = 0.75 \u00D7 \u221A(f\u2032c in psi) \u00D7 bw / fy"),
        eq("     = 0.75 \u00D7 0.07746 \u00D7 7.50 / 60"),
        eqResult("Av/s_min1 = 0.00726 in.\u00B2/in.    \u2190 GOVERNS"),
        blank(),
        P([B("Criterion 2 (Eq. 9.6.3.3b):")]),
        eq("Av/s = (50/1000) \u00D7 bw / fy = 0.050 \u00D7 7.50 / 60"),
        eqResult("Av/s_min2 = 0.00625 in.\u00B2/in."),
        blank(),
        P([B("Criterion 3 (Eq. 9.6.3.3c):")]),
        eq("Av/s = (Aps \u00D7 fpu)/(80 \u00D7 fy \u00D7 dp) \u00D7 \u221A(dp/bw)"),
        eq("     = (0.612 \u00D7 270)/(80 \u00D7 60 \u00D7 22.40) \u00D7 \u221A(22.40/7.50)"),
        eq("     = 165.24/107,520 \u00D7 1.728"),
        eqResult("Av/s_min3 = 0.00266 in.\u00B2/in."),
        blank(),
        P([B("Governing: "), R("Av/s_min = max(0.00726, 0.00625, 0.00266) = "), B("0.00726 in.\u00B2/in."), R(" (Criterion 1)")]),
        blank(),
        P([B("Stirrup spacing limits (ACI 318-19 Sec. 9.7.6.2.2):")]),
        eq("Basic:  s_max = min(3h/4, 24 in.) = min(21.0, 24) = 21.0 in."),
        eq("Tight trigger: Vs > 4\u221A(f\u2032c) \u00D7 bw \u00D7 dp = 52.06 kip"),
        eq("  Vs,req = 7.40 kip \u2264 52.06 kip  \u2192  basic limit applies"),

        // ══════════════════════════════════════════════════════════
        //  PART 5 \u2014 SUMMARY FOR ALL FOUR DESIGN SECTIONS
        // ══════════════════════════════════════════════════════════
        heading1("Part 5 \u2014 Summary at All Design Sections"),
        P([R("Values computed at four design sections. bw = 7.50 in., dp = 22.40 in. (0.8h governs throughout), \u03C6v = 0.75.")]),
        blank(),

        new Table({
          width: { size: 9648, type: WidthType.DXA },
          columnWidths: [960, 960, 960, 960, 1080, 1080, 960, 960, 960, 768],
          rows: [
            makeRow([
              makeCell("Section", 960, { shading: hdrShade, bold: true }),
              makeCell("x (ft)", 960, { shading: hdrShade, bold: true }),
              makeCell("Vu (kip)", 960, { shading: hdrShade, bold: true }),
              makeCell("Vci (kip)", 960, { shading: hdrShade, bold: true }),
              makeCell("Vcw (kip)", 1080, { shading: hdrShade, bold: true }),
              makeCell("Vc (kip)", 1080, { shading: hdrShade, bold: true }),
              makeCell("\u03C6Vc (kip)", 960, { shading: hdrShade, bold: true }),
              makeCell("Vs (kip)", 960, { shading: hdrShade, bold: true }),
              makeCell("Av/s", 960, { shading: hdrShade, bold: true }),
              makeCell("s (in.)", 768, { shading: hdrShade, bold: true }),
            ]),
            // Section 1: x = 18 in = 1.5 ft — Vci computed=102.36, capped to 65.07
            makeRow([
              makeCell("1", 960), makeCell("1.5", 960), makeCell("48.21", 960),
              makeCell("65.07\u2020", 960), makeCell("56.88", 1080), makeCell("56.88", 1080),
              makeCell("42.66", 960), makeCell("7.40", 960), makeCell("0.00726", 960),
              makeCell("21", 768),
            ]),
            // Section 2: x = 120 in = 10 ft — Mcr=0, Vci=max(24.47,22.12)=24.47
            makeRow([
              makeCell("2", 960), makeCell("10.0", 960), makeCell("34.78", 960),
              makeCell("24.47*", 960), makeCell("56.88", 1080), makeCell("24.47", 1080),
              makeCell("18.35", 960), makeCell("21.90", 960), makeCell("0.01629", 960),
              makeCell("13", 768),
            ]),
            // Section 3: x = 240 in = 20 ft — Mcr=0, floor governs: Vci=22.12
            makeRow([
              makeCell("3", 960), makeCell("20.0", 960), makeCell("18.97", 960),
              makeCell("22.12*", 960), makeCell("54.34", 1080), makeCell("22.12", 1080),
              makeCell("16.59", 960), makeCell("3.17", 960), makeCell("0.00726", 960),
              makeCell("21", 768),
            ]),
            // Section 4: x = 384 in = 32 ft (midspan) — Mcr=0, floor governs: Vci=22.12
            makeRow([
              makeCell("4", 960), makeCell("32.0", 960), makeCell("0.00", 960),
              makeCell("22.12*", 960), makeCell("54.34", 1080), makeCell("22.12", 1080),
              makeCell("16.59", 960), makeCell("0", 960), makeCell("0.00726", 960),
              makeCell("21", 768),
            ]),
          ],
        }),
        P([Ri("(*) Mcr = 0 at these sections; Vci governed by ACI floor (1.7\u03BB\u221Af\u2032c\u00D7bw\u00D7dp = 22.12 kip) or slightly above.", { size: 20 })]),
        P([Ri("(\u2020) Vci three-term = 102.36 kip, capped to Vci,max = 5.0\u03BB\u221Af\u2032c\u00D7bw\u00D7dp = 65.07 kip (CEE 530).", { size: 20 })]),
        blank(),
        P([B("Key observation: "), R("Vcw governs at Section 1 (near support). At Sections 2\u20134, Vci governs because Mcr = 0 drops the three-term Vci to the ACI floor (22.12 kip), which is well below Vcw. Section 2 requires significant stirrup reinforcement (Vs = 21.90 kip, s = 13 in.).")]),

        // ══════════════════════════════════════════════════════════
        //  PART 6 \u2014 DESIGN RECOMMENDATION
        // ══════════════════════════════════════════════════════════
        heading1("Part 6 \u2014 Design Recommendation"),
        P([R("Based on the analysis at all design sections:")]),
        blank(),
        P([B("Recommended stirrup design: #3 U-stirrups @ 13 in. in the support region (0\u201310 ft), transitioning to 21 in. elsewhere.")]),
        blank(),
        P([R("With the CEE 530 two-tier bounds (floor = 1.7, cap = 5.0), Vci governs at Sections 2\u20134 because Mcr = 0 drops the three-term Vci to the ACI floor (22.12 kip). Only at Section 1 does Vcw govern (56.88 kip). The critical section is Section 2 (x = 10 ft), where Vs = 21.90 kip requires Av/s = 0.01629 in.\u00B2/in., yielding s = 13 in. Beyond the drape point (x \u2265 20 ft), minimum Av/s governs and the spacing limit of 3h/4 = 21 in. controls.")]),

        // ══════════════════════════════════════════════════════════
        //  PART 7 \u2014 SHEAR DESIGN FIGURE
        // ══════════════════════════════════════════════════════════
        heading1("Part 7 \u2014 Shear Design Figure"),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { before: 120, after: 40 },
          children: [new ImageRun({
            type: "png", data: imgShear,
            transformation: { width: 580, height: 380 },
            altText: { title: "Shear Design", description: "4-panel shear design plots", name: "shear" },
          })],
        }),
        new Paragraph({
          alignment: AlignmentType.CENTER,
          spacing: { after: 120 },
          children: [Ri("Fig. 3 \u2014 (a) Shear demand vs. concrete capacity; (b) Required Vs; (c) Required Av/s; (d) Stirrup spacing along span.", { size: 20 })],
        }),

        // ══════════════════════════════════════════════════════════
        //  PART 8 \u2014 FINAL SUMMARY TABLE
        // ══════════════════════════════════════════════════════════
        heading1("Part 8 \u2014 Final Summary"),

        new Table({
          width: { size: 9648, type: WidthType.DXA },
          columnWidths: [3600, 2016, 2016, 2016],
          rows: [
            makeRow([
              makeCell("Check", 3600, { shading: hdrShade, bold: true }),
              makeCell("Demand", 2016, { shading: hdrShade, bold: true }),
              makeCell("Capacity", 2016, { shading: hdrShade, bold: true }),
              makeCell("Status", 2016, { shading: hdrShade, bold: true }),
            ]),
            makeRow([
              makeCell("dp \u2265 0.80h", 3600),
              makeCell("22.40 in.", 2016),
              makeCell("22.40 in.", 2016),
              makeCell("OK (=)", 2016, { bold: true }),
            ]),
            makeRow([
              makeCell("Vu at Section 1", 3600),
              makeCell("48.21 kip", 2016),
              makeCell("\u2014", 2016),
              makeCell("\u2014", 2016),
            ]),
            makeRow([
              makeCell("\u03C6Vc at Section 1", 3600),
              makeCell("\u2014", 2016),
              makeCell("42.66 kip", 2016),
              makeCell("\u2014", 2016),
            ]),
            makeRow([
              makeCell("Stirrups required?", 3600),
              makeCell("Vu > \u03C6Vc", 2016),
              makeCell("\u2014", 2016),
              makeCell("YES", 2016, { bold: true }),
            ]),
            makeRow([
              makeCell("Vs,req (Section 1)", 3600),
              makeCell("7.40 kip", 2016),
              makeCell("\u2014", 2016),
              makeCell("\u2014", 2016),
            ]),
            makeRow([
              makeCell("Vs \u2264 Vs,max (section size)", 3600),
              makeCell("7.40 kip", 2016),
              makeCell("104.1 kip", 2016),
              makeCell("PASS", 2016, { bold: true }),
            ]),
            makeRow([
              makeCell("\u03C6(Vc + Vs) \u2265 Vu", 3600),
              makeCell("48.21 kip", 2016),
              makeCell("53.22 kip", 2016),
              makeCell("PASS", 2016, { bold: true }),
            ]),
            makeRow([
              makeCell("Av/s_min (Crit. 1 governs)", 3600),
              makeCell("0.00726", 2016),
              makeCell("\u2014", 2016),
              makeCell("\u2014", 2016),
            ]),
            makeRow([
              makeCell("s at Section 1", 3600),
              makeCell("21 in.", 2016),
              makeCell("21 in.", 2016),
              makeCell("USE", 2016, { bold: true }),
            ]),
          ],
        }),

        blank(),
        heading2("ACI Code References"),
        new Table({
          width: { size: 9648, type: WidthType.DXA },
          columnWidths: [3600, 6048],
          rows: [
            makeRow([
              makeCell("ACI 318-19 Reference", 3600, { shading: hdrShade, bold: true }),
              makeCell("Description", 6048, { shading: hdrShade, bold: true }),
            ]),
            ...([
              ["Table 22.5.8.2", "Detailed method (Vci / Vcw)"],
              ["Eq. 22.5.8.3.1a", "Vci flexural-shear (floor: 1.7\u03BB\u221Af\u2032c\u00D7bw\u00D7dp; CEE 530 uses 5.0)"],
              ["Eq. 22.5.8.3.2a", "Mcr cracking moment"],
              ["Eq. 22.5.8.3.2", "Vcw web-shear"],
              ["Sec. 22.5.3.2", "Critical section at dp from support"],
              ["Table 9.6.3.3", "Minimum Av/s (3 criteria)"],
              ["Sec. 9.7.6.2.2", "Stirrup spacing limits (3h/4, 24 in.)"],
              ["Table 21.2.1", "\u03C6 = 0.75 for shear"],
              ["Sec. 5.3.1b", "Factored load combination (1.2D + 1.6L)"],
              ["Sec. 22.5.1.2", "Maximum Vs limit (8\u221Af\u2032c \u00D7 bw \u00D7 dp)"],
            ].map(([ref, desc]) => makeRow([
              makeCell(ref, 3600),
              makeCell(desc, 6048),
            ]))),
          ],
        }),
      ],
    },
  ],
});

// ── GENERATE ────────────────────────────────────────────────────────
Packer.toBuffer(doc).then(buffer => {
  const out = "ShearDesign_Report_v3.docx";
  fs.writeFileSync(out, buffer);
  console.log("Created:", out, "(" + (buffer.length/1024).toFixed(0) + " KB)");
});
