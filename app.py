from __future__ import annotations

import csv
import html
import io
import zipfile
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
import streamlit as st
import streamlit.components.v1 as components

from plasmid_verify.assembly import (
    TYPE_IIS_ENZYMES,
    AssemblyError,
    assemble_gibson,
    assemble_golden_gate,
)
from plasmid_verify.fasta import first_record, infer_target_id, to_fasta
from plasmid_verify.models import AssemblyResult, Mutation, VerificationResult
from plasmid_verify.verify import verify_consensus


st.set_page_config(
    page_title="CloningCompanion",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

st.markdown(
    """
<style>
  :root { --ink:#14261f; --muted:#607169; --paper:#fbfcf8; --green:#1c7854; --lime:#dff3df; --gold:#e7a93b; }
  .stApp { background: radial-gradient(circle at 85% 0%, #eaf4e9 0, transparent 25rem), #fbfcf8; }
  [data-testid="stSidebar"] { background:#12251e; }
  [data-testid="stSidebar"] * { color:#f5f8f3 !important; }
  [data-testid="stSidebar"] input { color:#14261f !important; }
  [data-testid="stSidebar"] [data-baseweb="input"] button svg { fill:#42564c !important; }
  .hero { padding:1.3rem 0 1.1rem; border-bottom:1px solid #ccd8d0; margin-bottom:1.1rem; }
  .eyebrow { color:#1c7854; font-weight:800; letter-spacing:.13em; text-transform:uppercase; font-size:.73rem; }
  .hero h1 { color:#14261f; font-size:2.45rem; line-height:1.05; margin:.3rem 0 .5rem; }
  .hero p { color:#607169; font-size:1.05rem; max-width:760px; margin:0; }
  .step { display:inline-flex; align-items:center; justify-content:center; width:1.65rem; height:1.65rem;
          border-radius:50%; background:#1c7854; color:white; font-weight:800; margin-right:.45rem; }
  .method-card { border:1px solid #cddbd2; border-radius:14px; padding:.9rem 1rem; background:rgba(255,255,255,.72); }
  .small-note { color:#607169; font-size:.84rem; }
  div[data-testid="stMetric"] { background:white; border:1px solid #d7e1da; padding:.8rem 1rem; border-radius:12px; }
  .verdict-pass { color:#166641; font-weight:850; }
  .verdict-review { color:#9b650d; font-weight:850; }
  .verdict-fail { color:#b43b38; font-weight:850; }
</style>
<div class="hero">
  <div class="eyebrow">Whole-plasmid sequencing · construct QC</div>
  <h1>CloningCompanion</h1>
  <p>Build the intended circular construct from a parent and ordered insert FASTAs, then explain every difference in the whole-plasmid consensus.</p>
</div>
""",
    unsafe_allow_html=True,
)


def _uploads_signature(files) -> Tuple[str, ...]:
    return tuple(file.name for file in (files or []))


def _default_insert_table(files) -> pd.DataFrame:
    rows = []
    order_by_target: Counter = Counter()
    for upload in files or []:
        target_id = infer_target_id(upload.name)
        order_by_target[target_id] += 1
        rows.append(
            {
                "insert_file": upload.name,
                "target_id": target_id,
                "assembly_order": order_by_target[target_id],
                "5prime_noncoding": 0,
                "3prime_noncoding": 0,
            }
        )
    return pd.DataFrame(rows)


def _default_consensus_table(files) -> pd.DataFrame:
    return pd.DataFrame(
        [{"consensus_file": upload.name, "target_id": infer_target_id(upload.name)} for upload in files or []]
    )


def _update_mapping_state(insert_files, consensus_files) -> None:
    insert_signature = _uploads_signature(insert_files)
    consensus_signature = _uploads_signature(consensus_files)
    if st.session_state.get("insert_signature") != insert_signature:
        st.session_state.insert_mapping = _default_insert_table(insert_files)
        st.session_state.insert_signature = insert_signature
    if st.session_state.get("consensus_signature") != consensus_signature:
        st.session_state.consensus_mapping = _default_consensus_table(consensus_files)
        st.session_state.consensus_signature = consensus_signature


def _mutation_rows(results: Dict[str, Tuple[AssemblyResult, VerificationResult]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for target_id, (_, result) in results.items():
        for mutation in result.mutations:
            row = {"target_id": target_id, "verdict": result.verdict}
            row.update(mutation.as_dict())
            rows.append(row)
    return rows


def _summary_rows(results: Dict[str, Tuple[AssemblyResult, VerificationResult]]) -> List[Dict[str, object]]:
    rows = []
    categories = ["backbone", "junction", "non-coding", "silent", "missense", "nonsense", "frameshift"]
    for target_id, (assembly, result) in results.items():
        counts = result.counts
        row = {
            "target_id": target_id,
            "verdict": result.verdict,
            "method": assembly.method,
            "target_bp": result.expected_length,
            "WPS_bp": result.observed_length,
            "mutations": len(result.mutations),
        }
        row.update({category: counts.get(category, 0) for category in categories})
        rows.append(row)
    return rows


REGION_COLORS = {
    "backbone": "#8da399",
    "junction": "#e7a93b",
    "non-coding": "#7c6db0",
    "coding": "#2a8c68",
}

CATEGORY_COLORS = {
    "backbone": "#657c71",
    "junction": "#d88a15",
    "non-coding": "#7c6db0",
    "silent": "#4a94b8",
    "missense": "#d15b45",
    "nonsense": "#a22f35",
    "stop loss": "#a22f35",
    "frameshift": "#711e45",
    "in-frame insertion": "#ca4f44",
    "in-frame deletion": "#ca4f44",
}


def _annotation_runs(assembly: AssemblyResult):
    if not assembly.annotations:
        return []
    runs = []
    start = 0
    current = assembly.annotations[0]
    for index, annotation in enumerate(assembly.annotations[1:], 1):
        if (annotation.region, annotation.source) != (current.region, current.source):
            runs.append((start, index, current))
            start, current = index, annotation
    runs.append((start, len(assembly.annotations), current))
    return runs


def _mutation_map(assembly: AssemblyResult, result: VerificationResult) -> str:
    width, height = 1100, 165
    left, usable = 40, width - 80
    length = max(1, len(assembly.sequence))
    segments = []
    for start, end, annotation in _annotation_runs(assembly):
        x = left + usable * start / length
        w = max(1.5, usable * (end - start) / length)
        label = html.escape(f"{annotation.source} · {annotation.region} · {start + 1}–{end}")
        color = REGION_COLORS.get(annotation.region, "#999")
        segments.append(
            f'<rect x="{x:.2f}" y="63" width="{w:.2f}" height="25" rx="3" fill="{color}"><title>{label}</title></rect>'
        )
    pins = []
    for mutation in result.mutations:
        x = left + usable * max(0, mutation.target_position - 1) / length
        color = CATEGORY_COLORS.get(mutation.category, "#b43b38")
        tip = html.escape(
            f"{mutation.target_position}: {mutation.category} · {mutation.target_bases}→{mutation.observed_bases} {mutation.protein_change}"
        )
        pins.append(
            f'<path d="M {x - 5:.2f} 53 L {x + 5:.2f} 53 L {x:.2f} 62 Z" fill="{color}" stroke="#fff" stroke-width="1"><title>{tip}</title></path>'
        )
        pins.append(f'<line x1="{x:.2f}" y1="43" x2="{x:.2f}" y2="53" stroke="{color}" stroke-width="2"/>')
    legend = []
    x = left
    for label, color in REGION_COLORS.items():
        legend.append(f'<rect x="{x}" y="118" width="11" height="11" rx="2" fill="{color}"/>')
        legend.append(f'<text x="{x + 16}" y="128" font-size="12" fill="#506159">{html.escape(label)}</text>')
        x += 120
    empty = "" if result.mutations else '<text x="550" y="38" text-anchor="middle" font-size="13" fill="#1c7854">No differences detected</text>'
    return f"""
    <div style="background:#fff;border:1px solid #d7e1da;border-radius:14px;padding:8px 10px;font-family:Inter,system-ui,sans-serif">
      <svg viewBox="0 0 {width} {height}" width="100%" role="img" aria-label="Linear target and mutation map">
        <text x="40" y="24" font-size="13" font-weight="700" fill="#14261f">TARGET MAP · {length:,} BP</text>
        {empty}{''.join(segments)}{''.join(pins)}
        <line x1="40" y1="94" x2="1060" y2="94" stroke="#cbd7d0"/>
        <text x="40" y="108" font-size="10" fill="#718078">1</text>
        <text x="1060" y="108" text-anchor="end" font-size="10" fill="#718078">{length:,}</text>
        {''.join(legend)}
      </svg>
    </div>"""


def _alignment_window(result: VerificationResult, mutation: Mutation, flank: int = 18) -> str:
    # Locate the alignment column corresponding to the 1-based target position.
    target_seen = 0
    center = 0
    for index, base in enumerate(result.aligned_target):
        if base != "-":
            target_seen += 1
        if target_seen >= mutation.target_position:
            center = index
            break
    start, end = max(0, center - flank), min(len(result.aligned_target), center + flank + 1)
    expected = result.aligned_target[start:end]
    observed = result.aligned_observed[start:end]
    marks = "".join("|" if a == b else "•" for a, b in zip(expected, observed))
    return f"Target    {expected}\n          {marks}\nObserved  {observed}"


with st.sidebar:
    st.markdown("### Run configuration")
    method = st.radio("Cloning method", ["Gibson", "Golden Gate"], horizontal=False)
    if method == "Gibson":
        arm_length = st.number_input("Terminal homology arm", min_value=10, max_value=60, value=15, step=1, help="Bases taken from each end of every insert FASTA.")
        enzyme_name = None
    else:
        enzyme_name = st.selectbox("Type IIS enzyme", list(TYPE_IIS_ENZYMES))
        enzyme = TYPE_IIS_ENZYMES[enzyme_name]
        st.caption(f"Recognition: {enzyme.recognition} · {enzyme.overhang_length}-nt fusion overhang")
        arm_length = 15
    st.markdown("---")
    st.caption("FASTA or SnapGene .dna · sequences are processed locally in this session")


setup_tab, results_tab, guide_tab = st.tabs(["1 · Inputs & assembly", "2 · Verification results", "Method notes"])

with setup_tab:
    st.markdown("### <span class='step'>1</span> Parent plasmid", unsafe_allow_html=True)
    parent_file = st.file_uploader(
        "Drop one circular parental plasmid sequence",
        type=["fa", "fasta", "fna", "fas", "dna"],
        accept_multiple_files=False,
        key="parent_upload",
    )

    st.markdown("### <span class='step'>2</span> Ordered inserts", unsafe_allow_html=True)
    st.caption("Name files with the target ID first, for example `T01_insert_01.dna`, `T01_insert_02.fasta`.")
    insert_files = st.file_uploader(
        "Drop insert sequences",
        type=["fa", "fasta", "fna", "fas", "dna"],
        accept_multiple_files=True,
        key="insert_uploads",
    )

    st.markdown("### <span class='step'>3</span> Whole-plasmid consensus", unsafe_allow_html=True)
    st.caption("Upload one assembled WPS consensus per target, such as `T01_consensus.dna`.")
    consensus_files = st.file_uploader(
        "Drop whole-plasmid sequencing consensuses",
        type=["fa", "fasta", "fna", "fas", "dna"],
        accept_multiple_files=True,
        key="consensus_uploads",
    )

    _update_mapping_state(insert_files, consensus_files)
    if insert_files:
        st.markdown("#### Check insert pairing and order")
        st.caption("Non-coding values mark bases between each terminal junction and the coding frame.")
        st.session_state.insert_mapping = st.data_editor(
            st.session_state.insert_mapping,
            hide_index=True,
            use_container_width=True,
            disabled=["insert_file"],
            column_config={
                "target_id": st.column_config.TextColumn(required=True),
                "assembly_order": st.column_config.NumberColumn(min_value=1, step=1, required=True),
                "5prime_noncoding": st.column_config.NumberColumn("5′ non-coding (bp)", min_value=0, step=1),
                "3prime_noncoding": st.column_config.NumberColumn("3′ non-coding (bp)", min_value=0, step=1),
            },
            key="insert_mapping_editor",
        )
    if consensus_files:
        st.markdown("#### Check WPS-to-target pairing")
        st.session_state.consensus_mapping = st.data_editor(
            st.session_state.consensus_mapping,
            hide_index=True,
            use_container_width=True,
            disabled=["consensus_file"],
            column_config={"target_id": st.column_config.TextColumn(required=True)},
            key="consensus_mapping_editor",
        )

    ready = bool(parent_file and insert_files and consensus_files)
    run = st.button("Build targets & verify WPS", type="primary", use_container_width=True, disabled=not ready)
    if run:
        errors: List[str] = []
        results: Dict[str, Tuple[AssemblyResult, VerificationResult]] = {}
        try:
            parent = first_record(parent_file).sequence
        except Exception as exc:
            parent = ""
            errors.append(f"Parent sequence: {exc}")

        insert_by_name = {upload.name: upload for upload in insert_files}
        consensus_by_name = {upload.name: upload for upload in consensus_files}
        insert_rows = st.session_state.insert_mapping.to_dict("records")
        consensus_rows = st.session_state.consensus_mapping.to_dict("records")
        progress = st.progress(0, text="Preparing targets…")
        for number, consensus_row in enumerate(consensus_rows, 1):
            target_id = str(consensus_row["target_id"]).strip()
            try:
                matching = [row for row in insert_rows if str(row["target_id"]).strip() == target_id]
                if not matching:
                    raise AssemblyError("No inserts are mapped to this target")
                matching.sort(key=lambda row: int(row["assembly_order"]))
                inserts = []
                noncoding = {}
                for row in matching:
                    upload = insert_by_name[row["insert_file"]]
                    record = first_record(upload)
                    inserts.append((row["insert_file"], record.sequence))
                    noncoding[row["insert_file"]] = (
                        int(row.get("5prime_noncoding", 0) or 0),
                        int(row.get("3prime_noncoding", 0) or 0),
                    )
                if method == "Gibson":
                    assembly = assemble_gibson(parent, inserts, target_id, int(arm_length), noncoding)
                else:
                    assembly = assemble_golden_gate(parent, inserts, target_id, enzyme_name, noncoding)
                observed = first_record(consensus_by_name[consensus_row["consensus_file"]]).sequence
                verification = verify_consensus(assembly, observed)
                results[target_id] = (assembly, verification)
            except Exception as exc:
                errors.append(f"{target_id or consensus_row['consensus_file']}: {exc}")
            progress.progress(number / len(consensus_rows), text=f"Processed {number} of {len(consensus_rows)} targets")
        progress.empty()
        st.session_state.verification_results = results
        st.session_state.run_errors = errors
        if results:
            st.success(f"Verified {len(results)} target{'s' if len(results) != 1 else ''}. Open the results tab.")
        if errors:
            st.error("Some targets could not be processed:\n\n" + "\n\n".join(f"• {error}" for error in errors))

with results_tab:
    results = st.session_state.get("verification_results", {})
    errors = st.session_state.get("run_errors", [])
    if not results:
        st.info("Upload the parent, inserts, and WPS consensus files, then run verification.")
    else:
        summary = pd.DataFrame(_summary_rows(results))
        verdict_counts = Counter(row["verdict"] for row in summary.to_dict("records"))
        metric_cols = st.columns(4)
        metric_cols[0].metric("Targets", len(summary))
        metric_cols[1].metric("Pass", verdict_counts.get("PASS", 0))
        metric_cols[2].metric("Review", verdict_counts.get("REVIEW", 0))
        metric_cols[3].metric("Fail", verdict_counts.get("FAIL", 0))
        st.dataframe(summary, hide_index=True, use_container_width=True)

        download_cols = st.columns(3)
        download_cols[0].download_button(
            "Download summary CSV",
            summary.to_csv(index=False).encode(),
            "plasmid_verification_summary.csv",
            "text/csv",
            use_container_width=True,
        )
        mutation_frame = pd.DataFrame(_mutation_rows(results))
        download_cols[1].download_button(
            "Download mutations CSV",
            mutation_frame.to_csv(index=False).encode(),
            "plasmid_verification_mutations.csv",
            "text/csv",
            use_container_width=True,
        )
        target_zip = io.BytesIO()
        with zipfile.ZipFile(target_zip, "w", zipfile.ZIP_DEFLATED) as archive:
            for target_id, (assembly, _) in results.items():
                archive.writestr(f"{target_id}_target.fasta", to_fasta(target_id, assembly.sequence))
        download_cols[2].download_button(
            "Download target FASTAs",
            target_zip.getvalue(),
            "generated_targets.zip",
            "application/zip",
            use_container_width=True,
        )

        st.markdown("### Target detail")
        selected = st.selectbox("Target", list(results))
        assembly, verification = results[selected]
        verdict_class = verification.verdict.lower()
        st.markdown(
            f"<span class='verdict-{verdict_class}'>{verification.verdict}</span> · "
            f"{html.escape(assembly.method)} · {len(verification.mutations)} difference(s) · "
            f"{verification.orientation}, rotated {verification.rotation:,} bp",
            unsafe_allow_html=True,
        )
        components.html(_mutation_map(assembly, verification), height=188, scrolling=False)

        with st.expander("How this target was assembled"):
            for note in assembly.notes:
                st.write(f"• {note}")

        if not verification.mutations:
            st.success("The WPS consensus matches the generated target across the full circular plasmid.")
        else:
            mutation_rows = [mutation.as_dict() for mutation in verification.mutations]
            st.dataframe(pd.DataFrame(mutation_rows), hide_index=True, use_container_width=True)
            st.markdown("#### Sequence-level evidence")
            for mutation in verification.mutations:
                title = (
                    f"{mutation.target_position:,} · {mutation.category} · "
                    f"{mutation.target_bases} → {mutation.observed_bases}"
                )
                if mutation.protein_change:
                    title += f" · {mutation.protein_change}"
                with st.expander(title):
                    st.code(_alignment_window(verification, mutation), language="text")
                    st.caption(f"Region: {mutation.region} · source: {mutation.source}")
        if errors:
            with st.expander(f"{len(errors)} target(s) with input or assembly errors"):
                for error in errors:
                    st.error(error)

with guide_tab:
    left, right = st.columns(2)
    with left:
        st.markdown("### Gibson assembly")
        st.markdown(
            "The first and last terminal bases of each insert are treated as homology arms. "
            "The app searches the circular parent in both orientations, replaces the shortest interval "
            "bracketed by those arms, and repeats in assembly order for multi-insert targets."
        )
        st.code(
            ">T01_insert_01\n"
            "[15-nt LEFT ARM][insert payload][15-nt RIGHT ARM]",
            language="text",
        )
    with right:
        st.markdown("### Golden Gate assembly")
        st.markdown(
            "Each insert must include an inward-facing pair of the selected Type IIS sites. "
            "The parent must contain a matching destination cassette. Recognition sites and spacers "
            "are removed; ordered parts are accepted only when every fusion overhang is compatible."
        )
        st.code(
            "5′ [Type IIS →][fusion][payload][fusion][← Type IIS] 3′",
            language="text",
        )
    st.markdown("### Verdict rules")
    st.markdown(
        "**PASS** means no sequence differences. **REVIEW** covers only silent, non-coding, or backbone changes. "
        "**FAIL** marks junction changes or protein-altering coding changes, including missense, nonsense, "
        "stop-loss, frameshift, and in-frame indels. Ambiguous `N` bases are retained and reported as differences."
    )
