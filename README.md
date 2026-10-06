# CloningCompanion

CloningCompanion generates intended circular plasmids from a parental plasmid and ordered insert sequences, then compares each generated target with one whole-plasmid sequencing (WPS) consensus. Its native macOS app also designs Golden Gate PCR and synthesis fragments directly from an annotated SnapGene map. Both apps accept FASTA and SnapGene `.dna` inputs.

**[Download for Apple silicon Macs](https://github.com/benno-orr/cloning-companion/releases/latest)** · macOS 12+

This repository contains application source and downloadable releases, not user sequence files or projects.
Do not post confidential sequences or project data in public issues. This is research software; review generated designs before ordering or experimental use.

## Run

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Native macOS application

The project also includes a native macOS shell that embeds Python and Biopython in a WebKit-based `.app`. It does not launch a browser or run a localhost server. Sequence and project selection use macOS open panels; CSV, ZIP, and project exports use macOS save panels.

The native application includes:

- File and Help menus
- About and workflow help screens
- `.plasmidverify` project documents
- Recent-project history in Application Support
- Native FASTA/SnapGene selection and report export
- The same tested Gibson, Golden Gate, circular alignment, and mutation-classification engine

Build locally:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements-macos.txt
./scripts/build_macos.sh
```

Outputs:

```text
dist/CloningCompanion.app
dist/CloningCompanion-1.8.0.dmg
```

When `SIGN_IDENTITY` is omitted, the build script uses an ad-hoc signature suitable for local testing. For public distribution, use a Developer ID Application certificate:

```bash
SIGN_IDENTITY="Developer ID Application: Your Name (TEAMID)" ./scripts/build_macos.sh
```

### In-app local updates

Version 1.3.0 adds **Check for Updates → Update & Restart** in the sidebar and Help menu.
Each successful build publishes its version and bundle checksum to this Mac’s Application Support folder.
Use `BUILD_DMG=0 ./scripts/build_macos.sh` to build and publish without creating a DMG.
The updater stages and validates the app (checksum, bundle identity, code signature), saves unsaved project inputs,
waits for the running app to exit, replaces it, and reopens it. A `CloningCompanion.previous-*.app` backup remains
beside the installed app, or under Application Support/CloningCompanion/Updates when the app’s parent folder is protected.
For a writable app in a protected parent folder, the updater swaps its Contents on the same filesystem.
If replacement or launching fails, the old bundle is restored.
Inputs restore automatically. Completed runs created in 1.7.0 or later remain in the local app library across restarts.
When upgrading from 1.6.0, download any temporary results you need before restarting, or regenerate them afterward.
Previously downloaded files remain unchanged; verification results can be rerun.
Recovery files and helper logs are under `~/Library/Application Support/CloningCompanion/Updates/install-*`.
Pre-1.3.0 installations need a one-time app replacement to acquire the updater; subsequent updates need no DMG.
Read-only/translocated apps must first be moved to Applications.

### GitHub updates and downloads

[Download the latest release](https://github.com/benno-orr/cloning-companion/releases/latest).
The desktop build currently supports **Apple silicon, macOS 12 or later**.
Unzip the download and move CloningCompanion.app to Applications before launching.
The current builds are ad-hoc signed, not Developer ID signed or Apple-notarized; macOS may require explicit
approval on first launch. Ed25519 update signatures are separate from Apple signing/notarization.

Version 1.4.0 and later checks this public repository’s GitHub Releases with no GitHub account required.
It verifies `latest.json` with an embedded Ed25519 public key, then verifies the downloaded ZIP’s signed hash,
size, extracted bundle checksum and macOS code signature before staging it. No sequence data is uploaded.
Newer locally built versions remain available on the developer’s Mac, taking priority over GitHub releases.

To prepare a release: run `BUILD_DMG=0 ./scripts/build_macos.sh`, then
`.venv/bin/python scripts/prepare_github_release.py dist/CloningCompanion.app`.
Upload only the versioned app ZIP and `latest.json` from `dist/github-release-VERSION/` as assets on release `vVERSION`.
Create a draft, upload both assets, then publish it. Do not overwrite published release assets.
The private signing key stays outside the repository in Application Support/CloningCompanion/Signing;
back it up securely. Only the public key in `mac_app/update_identity.py` belongs in source control.
Do not rotate that key without a planned trust-key transition for existing installations.

Version 1.5.1 restores translation tracks for pasted replacement variants when the input
explicitly annotates the entire core as coding. Translations are recalculated from the
selected DNA using the template reading frame and color, and labelled with the variant
name. Internal template domains and mutation labels are not inferred on replacements.
Linear-map callouts have individually fitted connectors and no redundant junction headings.

Store notarization credentials once using Apple's interactive prompt, then submit and staple the DMG:

```bash
xcrun notarytool store-credentials PlasmidVerifyNotary
NOTARY_PROFILE=PlasmidVerifyNotary ./scripts/notarize_macos.sh
```

The current PyInstaller target is Apple silicon (`arm64`). Producing a universal Intel/Apple-silicon build requires a universal Python distribution and universal wheels for compiled dependencies.

Inputs can be FASTA (`.fa`, `.fasta`, `.fna`, or `.fas`) or SnapGene (`.dna`). Formats may be mixed within a project. Use a shared target prefix:

```text
T01_consensus.fasta
T01_insert_01.dna
T01_insert_02.fasta
T02_consensus.fasta
T02_insert_01.fasta
```

The inferred pairing and multi-insert order are editable before the run. SnapGene files are read for their sequence; feature annotations are not imported. Generated targets still export as FASTA.

## Assembly models and design

### Unique assembled SnapGene products

The design workflow exports every combination of variable-fragment options (up to 10,000 combinations per run).
Each distinct double-stranded circular DNA sequence gets one annotated `.dna` file in the reported plasmid folder.
Equivalent rotations, reverse complements, and identically sequenced aliases share a file. `plasmids.tsv` lists
unique products, and `plasmid_combinations.tsv` maps every option combination to its file. The app lists the files
under **Assembled plasmids**, with individual **Show file** buttons. A rerun preserves previous plasmid folders.
Order fragments remain one per named option, not one per plasmid. Schematic graphics depict the first combination.

Native feature colors and top/bottom DNA-strand color ranges take precedence over the default HSV palette.
Unchanged annotations are projected into each product; changed variant cores inherit their parent fragment’s colors
without inheriting stale protein translations. New variant DNA uses the dominant parent-core color independently
for each strand; retained linker/flank bases keep their original colors. HSV (30° × n, 50%, 100%), backbone n=0,
is used only where no input color is available. Every product is read back and checked after writing.

- **Gibson:** the first and last 15 nt of each insert (configurable) are located in the circular parent in either orientation. The shortest bracketed parental interval is replaced by the full insert. Multi-insert assemblies are applied in the displayed order.
- **Golden Gate:** supports BsaI, BsmBI/Esp3I, BbsI, SapI, and AarI. The app finds an inward-facing site pair in the parent and each insert, removes recognition sites and spacers, checks all fusion overhangs, and produces the scarless target.

Golden Gate insert sequences must include the Type IIS recognition sites and fusion overhangs. The parent sequence must include the destination cassette. Internal occurrences of the selected enzyme can make a design ambiguous; the app reports incompatible or missing site arrangements rather than guessing.

The native app's **Golden Gate design** workspace takes one annotated SnapGene `.dna` map—no YAML needed. Label a fixed fragment core `[name]` and a variable fragment core `{name}`. The bracketed core overlapping `pUC ori` or `ori pUC` is automatically the backbone. Their order on the circular map is the assembly order. Add exactly one `(name)` or `-` feature touching each adjacent core boundary; this marks where that junction may be placed. Each highlighted junction region must be at least 4 bp; the app selects compatible non-palindromic 4 bp fusions inside those regions. Other ordinary annotations are ignored. The detected backbone is output as PCR and other cores as synthesis fragments. It writes `junctions.tsv`, `candidate_options.tsv`, `fragments.tsv`, `assembly_recipe.tsv`, `pcr_primers.tsv`, and `synthesis_order.fasta`.

The output now also includes `assembly_schematic.dna`: an editable planning map with gray `[fixed]` / `{variable}` blocks, separate top/bottom strand colors, staggered fusion coloring, and outward-facing Type IIS sites separated by gray N gaps. Original subfeatures are projected onto unchanged scaffold-derived cores; replaced variant cores do not inherit potentially stale subfeatures. The N gaps and duplicate fusions are schematic only and are excluded from all order DNA. `assembled_product_map.dna` is the separate scarless representative product; it contains no visual gaps or adapters.

Each junction also has a high-resolution `junction_screencaps/*_clean.png` figure: fragment ends face each other across white space, retaining strand colors, gray paired overhang bases, fragment bars, and available protein annotations. Enzyme adapters and spacer bases are hidden in this view. `junction_clean_overview.png` combines all junctions, starting with the backbone-to-insert junction. The detailed enzyme views and final assembled-junction images are also exported.

The **Joined junctions** gallery shows complementary ends aligned across a stepped cut line. `junction_screencaps/*_joined.png` and `junction_joined_overview.png` show each four-base fusion once, with both strands colored by native core membership, including fusions that cross a core boundary. Copies appended to a neighboring fragment retain their source core's color. Connecting DNA outside labelled cores is assigned to the downstream core for coloring. The same provenance is retained for pasted variants of different lengths. Noncoding bases are lowercase when native translated-feature annotations are available.

Drop an annotated `.dna` file onto the design-map input or click to choose one, name the project, then run. Starting in 1.8.0, completed designs are saved automatically **as projects inside the app**, including the exact input SnapGene map, enzyme, pasted variant lists, graphics and output bytes. **Graphics & orders → Saved projects** restores the inputs and outputs together. The embedded map can be used for reruns even after the original file is moved, deleted or edited. New successful runs create revisions of the same project; the picker shows its newest revision, while prior revisions remain in storage. Use New Project for a separate project. Storage stays local in `~/Library/Application Support/CloningCompanion/DesignLibrary/designs.sqlite3`; renderer staging is cleaned up and no output folder is exported automatically. Older 1.7 output-only runs remain accessible but cannot restore inputs that were never stored. Nothing is uploaded, and existing exported files are untouched.

Use **Download outputs** to export individual SnapGene plasmids, synthesis TSV/FASTA, PCR primers, schematics and reports, or ZIP bundles. **Save Project** optionally creates a portable version-2 `.plasmidverify` document with the Golden Gate input map and matching generated outputs embedded, not merely links to the current Mac's library. Opening it restores previews and downloads without regeneration and supports rerunning on another Mac. Saving rejects mismatched inputs and outputs, uses atomic file replacement, and leaves the current project association unchanged if Save As is cancelled. Legacy version-1 project files still open; verification input files remain external references as before. The CLI retains its explicit output-directory behavior. Version 1.6.0's memory-only results must be downloaded before upgrading or regenerated to enter the library.

For each `{variable core}`, paste a list of named core DNA sequences in its variant box: two tab-separated columns, CSV (optional `name` and `dna_sequence` headers), or FASTA. The app adds the selected shared fusions and enzyme adapters to every version. Fixed pieces are output once. `synthesis_order.tsv` and `synthesis_order.fasta` contain one synthesis fragment per version; PCR variants get individual primer pairs in `pcr_primers.tsv`. Choose one option per assembly position in `assembly_recipe.tsv`. The SnapGene assembled map and junction PNGs represent the first option from each list, not all combinations. Save Project preserves the pasted lists. Empty lists use the core from the original map, which must contain unambiguous DNA. These checks do not predict experimental ligation fidelity or substitute for vendor/manufacturability review.

The **Linear map** view displays a to-scale linearized plasmid with alternating upper/lower junction blowups, downloadable as `plasmid_with_junction_blowups.png`. **Fragment ends** shows the separate pre-digestion ends. Click an image to enlarge it or view it at actual size; each has a PNG download button.

Regions labelled `-` are fixed scaffold-derived sequence, never user-supplied library elements. Paste variable cores without any `-` bases, including marked portions overlapping a core label. The designer restores these bases automatically and distributes them between the two flanking fragments around the chosen four-base overlap. `automatic_junction_bases.tsv` and the in-app report show the full region and each fragment-end contribution. The assembled region contains every native base once. Older pasted lists containing overlapping `-` bases must have those marked bases removed before reuse.

## Mutation interpretation

Annotated-map design now cuts at native fusion coordinates rather than adding four bases to labelled cores. A junction may split a linker at any allowed position; all other native linker/flanking bases are retained. Without variants, the assembled product must exactly match the uploaded circular plasmid. With pasted variants, full labelled cores are replaced and all remaining native DNA is preserved. This independent sequence check runs before writing order files. Explicit YAML plans still declare fusion ownership; a repeated end fusion is rejected unless explicitly authorized with `allow_repeated_fusion`. Regenerate designs exported with versions 1.1.12 and earlier: they could duplicate owned fusion bases and omit junction-region/linker bases.

WPS sequences are automatically reverse-complemented and rotated to the generated circular target. A full-length global alignment reports substitutions and indels. Per-base target provenance divides events into backbone, junction, non-coding, and coding regions. Coding substitutions are translated as silent, missense, nonsense, or stop-loss; coding indels are classified as frameshift or in-frame.

Verdicts:

- `PASS`: exact full-plasmid match.
- `REVIEW`: only silent, non-coding, or backbone differences.
- `FAIL`: protein-altering coding or junction differences.

Results include a visual mutation map, alignment context for each event, summary and mutation CSVs, and generated target FASTAs.

## Test

```bash
pytest -q
```
