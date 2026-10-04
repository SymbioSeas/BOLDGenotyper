"""
Quality Control and Analysis Reports

This module generates summary reports and QC metrics for the haplotyping pipeline,
including taxonomy conflict reports, assignment summaries, and haplotype
characterizations.
"""

import logging
import json
from html import escape as html_escape
from pathlib import Path
from typing import Dict, Optional
import pandas as pd
import numpy as np

logger = logging.getLogger(__name__)


SUMMARY_COLUMNS = [
    "n_input", "n_unique", "n_no_sequence", "n_failed_qc", "n_qc_pass",
    "n_excluded_core", "n_core_ok", "n_assigned", "n_assigned_corrected",
    "n_excluded_error", "n_excluded_error_tied", "n_excluded_error_far", "n_haplotypes",
    "pct_assigned", "pct_assigned_exact",
    "n_species_labels_input", "n_species_labels_species_level", "n_species_labels_qc_passed",
    "n_species_haplotype_assigned", "n_species_confident",
]
SPECIES_COUNT_COLUMNS = [c for c in SUMMARY_COLUMNS if c.startswith("n_species_")]

# Record statuses (record_status.STATUSES) and their summary columns
_STATUS_COLUMNS = {
    "assigned": "n_assigned",
    "assigned_corrected": "n_assigned_corrected",
    "excluded_error_singleton": "n_excluded_error",
    "excluded_core_coverage": "n_excluded_core",
    "failed_qc": "n_failed_qc",
    "no_sequence": "n_no_sequence",
}


def compute_assignment_summary(status_df, n_input, n_haplotypes, species_counts,
                               n_excluded_error_tied=0, n_excluded_error_far=0):
    """Assignment summary with explicit denominators.

    Percentages are of QC-passed records (n_qc_pass = unique records minus
    no_sequence and failed_qc). Returns a dict with SUMMARY_COLUMNS keys.
    """
    counts = status_df["status"].value_counts()

    def c(status):
        return int(counts.get(status, 0))

    n_unique = int(len(status_df))
    n_qc_pass = n_unique - c("no_sequence") - c("failed_qc")
    summary = {
        "n_input": int(n_input if n_input is not None else n_unique),
        "n_unique": n_unique,
        "n_no_sequence": c("no_sequence"),
        "n_failed_qc": c("failed_qc"),
        "n_qc_pass": n_qc_pass,
        "n_excluded_core": c("excluded_core_coverage"),
        "n_core_ok": n_qc_pass - c("excluded_core_coverage"),
        "n_assigned": c("assigned"),
        "n_assigned_corrected": c("assigned_corrected"),
        "n_excluded_error": c("excluded_error_singleton"),
        "n_excluded_error_tied": int(n_excluded_error_tied),
        "n_excluded_error_far": int(n_excluded_error_far),
        "n_haplotypes": int(n_haplotypes),
        "pct_assigned": round((c("assigned") + c("assigned_corrected")) / n_qc_pass * 100, 2) if n_qc_pass else 0.0,
        "pct_assigned_exact": round(c("assigned") / n_qc_pass * 100, 2) if n_qc_pass else 0.0,
    }
    for key in SPECIES_COUNT_COLUMNS:
        summary[key] = int((species_counts or {}).get(key, 0))
    return {k: summary[k] for k in SUMMARY_COLUMNS}


def generate_assignment_summary(
    annotated_csv: str,
    diagnostics_csv: str,
    output_csv: str,
    n_input: Optional[int] = None,
    n_haplotypes: Optional[int] = None,
    species_counts: Optional[Dict[str, int]] = None,
    n_excluded_error_tied: int = 0,
    n_excluded_error_far: int = 0,
) -> Dict[str, object]:
    """
    Write the assignment summary: one row with explicit denominators.

    Parameters
    ----------
    annotated_csv : str
        Path to the annotated CSV (kept for compatibility; statuses come from
        the diagnostics table)
    diagnostics_csv : str
        Path to the diagnostics CSV: one row per unique input record with
        'status' and 'assigned_haplotype'
    output_csv : str
        Path for output summary CSV
    n_input : int, optional
        Data rows in the input TSV, before duplicate removal (default: the
        number of unique records)
    n_haplotypes : int, optional
        Number of retained haplotypes (default: distinct assigned haplotypes)
    species_counts : dict, optional
        taxon_labels.species_count_summary() output
    n_excluded_error_tied, n_excluded_error_far : int, optional
        Error-singleton records left unassigned because they are equidistant
        from two or more haplotypes, or because no retained haplotype is within
        the error-filter threshold

    Returns
    -------
    dict
        The summary (SUMMARY_COLUMNS keys)
    """
    out = Path(output_csv)
    out.parent.mkdir(parents=True, exist_ok=True)

    diag = pd.read_csv(diagnostics_csv)
    if n_haplotypes is None:
        n_haplotypes = diag["assigned_haplotype"].dropna().nunique() if "assigned_haplotype" in diag else 0

    summary = compute_assignment_summary(diag, n_input, n_haplotypes, species_counts,
                                         n_excluded_error_tied, n_excluded_error_far)
    pd.DataFrame([summary], columns=SUMMARY_COLUMNS).to_csv(out, index=False)
    logger.info(f"Assignment summary saved: {out}")
    return summary


# ============================================================================
# HTML Report Generation
# ============================================================================

# HTML Report CSS
HTML_REPORT_CSS = """
<style>
    :root {
        --primary-color: #2c3e50;
        --secondary-color: #3498db;
        --success-color: #27ae60;
        --warning-color: #f39c12;
        --danger-color: #e74c3c;
        --light-bg: #ecf0f1;
        --border-color: #bdc3c7;
    }

    * {
        margin: 0;
        padding: 0;
        box-sizing: border-box;
    }

    body {
        font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
        line-height: 1.6;
        color: #333;
        background-color: #f5f5f5;
        margin: 0;
        padding: 0;
    }

    .container {
        display: flex;
        min-height: 100vh;
        background-color: white;
    }
    
    .geo-map {
        position: relative;
        width: 100%;
        max-width: 100%;
        height: 500px;
        border: 1px solid var(--border-color);
        border-radius: 6px;
        background-color: #f8f9fa;
        margin-top: 15px;
        margin-bottom: 25px;
        overflow: hidden;
    }

    .geo-map-svg {
        width: 100%;
        height: 100%;
    }

    .geo-point {
        fill: var(--secondary-color);
        stroke: white;
        stroke-width: 1;
        opacity: 0.8;
        cursor: pointer;
    }

    .geo-point:hover {
        opacity: 1.0;
    }

    .geo-tooltip {
        position: absolute;
        display: none;
        pointer-events: none;
        background-color: rgba(0, 0, 0, 0.75);
        color: white;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.8em;
        z-index: 10;
        white-space: nowrap;
    }
    
    .header {
        background: linear-gradient(135deg, var(--primary-color), var(--secondary-color));
        color: white;
        padding: 30px 40px;
        text-align: center;
    }

    .header h1 {
        font-size: 2em;
        margin-bottom: 8px;
        font-style: italic;
    }

    .header .subtitle {
        font-size: 1em;
        opacity: 0.9;
        margin-bottom: 5px;
    }

    .header .timestamp {
        font-size: 0.85em;
        opacity: 0.8;
    }

    .quick-stats {
        background-color: var(--light-bg);
        padding: 20px 40px;
        display: flex;
        justify-content: space-around;
        flex-wrap: wrap;
        border-bottom: 3px solid var(--secondary-color);
    }

    .quick-stat {
        text-align: center;
        padding: 10px 20px;
    }

    .quick-stat .number {
        font-size: 2em;
        font-weight: bold;
        color: var(--secondary-color);
    }

    .quick-stat .label {
        font-size: 0.9em;
        color: #666;
        margin-top: 5px;
    }

    .content {
        padding: 40px;
        flex: 1;
    }

    .section {
        margin-bottom: 50px;
    }

    .section h2 {
        color: var(--primary-color);
        border-bottom: 2px solid var(--secondary-color);
        padding-bottom: 10px;
        margin-bottom: 20px;
        font-size: 1.8em;
    }

    .section h3 {
        color: var(--primary-color);
        margin-top: 30px;
        margin-bottom: 15px;
        font-size: 1.3em;
    }

    /* Left Sidebar Navigation */
    .sidebar {
        width: 300px;
        background: linear-gradient(180deg, var(--primary-color), #1a252f);
        color: white;
        position: fixed;
        height: 100vh;
        overflow-y: auto;
        padding: 20px 0;
        box-shadow: 2px 0 10px rgba(0,0,0,0.1);
        z-index: 1000;
    }

    .sidebar-header {
        padding: 0 20px 20px 20px;
        border-bottom: 1px solid rgba(255,255,255,0.2);
        margin-bottom: 20px;
    }

    .sidebar-title {
        font-size: 1.1em;
        font-weight: 600;
        margin-bottom: 5px;
        font-style: italic;
    }

    .sidebar-version {
        font-size: 0.75em;
        opacity: 0.7;
    }

    .sidebar-nav {
        list-style: none;
        padding: 0;
        margin: 0;
    }

    .sidebar-nav-item {
        margin: 0;
    }

    .sidebar-nav-link {
        display: block;
        padding: 12px 20px;
        text-decoration: none;
        color: rgba(255,255,255,0.9);
        border-left: 3px solid transparent;
        transition: all 0.2s ease;
        font-size: 0.95em;
    }

    .sidebar-nav-link:hover {
        background-color: rgba(255,255,255,0.1);
        border-left-color: var(--secondary-color);
        color: white;
    }

    .sidebar-nav-link.active {
        background-color: rgba(255,255,255,0.15);
        border-left-color: var(--secondary-color);
        color: white;
        font-weight: 600;
    }

    /* Main Content Area */
    .main-wrapper {
        flex: 1;
        margin-left: 300px;
        display: flex;
        flex-direction: column;
    }
    
    .metric-grid {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(310px, 1fr));
        gap: 20px;
        margin-bottom: 30px;
    }

    .metric-card {
        background-color: white;
        border: 1px solid var(--border-color);
        border-radius: 6px;
        padding: 20px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
    }

    .metric-card .title {
        font-size: 0.85em;
        color: #666;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        margin-bottom: 10px;
    }

    .metric-card .value {
        font-size: 2em;
        font-weight: bold;
        color: var(--secondary-color);
    }

    .metric-card .description {
        font-size: 0.9em;
        color: #888;
        margin-top: 5px;
    }

    .table-container {
        overflow-x: auto;
        margin: 20px 0;
    }

    table {
        width: 100%;
        border-collapse: collapse;
        font-size: 0.9em;
    }

    thead {
        background-color: var(--primary-color);
        color: white;
    }

    th {
        padding: 12px;
        text-align: left;
        font-weight: 600;
    }

    td {
        padding: 10px 12px;
        border-bottom: 1px solid #eee;
    }

    tbody tr:hover {
        background-color: #f8f9fa;
    }

    tbody tr:nth-child(even) {
        background-color: #fafafa;
    }

    .alert {
        padding: 15px;
        border-radius: 6px;
        margin: 20px 0;
        border-left: 4px solid;
    }

    .alert-info {
        background-color: #d1ecf1;
        border-color: #17a2b8;
        color: #0c5460;
    }

    .alert-warning {
        background-color: #fff3cd;
        border-color: #ffc107;
        color: #856404;
    }

    .alert-success {
        background-color: #d4edda;
        border-color: #28a745;
        color: #155724;
    }

    .badge {
        display: inline-block;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.85em;
        font-weight: 600;
    }

    .badge-success {
        background-color: var(--success-color);
        color: white;
    }

    .badge-warning {
        background-color: var(--warning-color);
        color: white;
    }

    .badge-danger {
        background-color: var(--danger-color);
        color: white;
    }

    .badge-info {
        background-color: var(--secondary-color);
        color: white;
    }

    .footer {
        background-color: var(--light-bg);
        padding: 30px 40px;
        text-align: center;
        border-top: 1px solid var(--border-color);
        font-size: 0.9em;
        color: #666;
    }

    .footer a {
        color: var(--secondary-color);
        text-decoration: none;
    }

    .footer a:hover {
        text-decoration: underline;
    }

    @media print {
        body {
            background-color: white;
            padding: 0;
        }
        .container {
            box-shadow: none;
        }
        .section {
            page-break-inside: avoid;
        }
    }

    @media (max-width: 768px) {
        .sidebar {
            width: 200px;
        }
        .main-wrapper {
            margin-left: 200px;
        }
        .header h1 {
            font-size: 1.5em;
        }
        .quick-stats {
            flex-direction: column;
            padding: 15px 20px;
        }
        .metric-grid {
            grid-template-columns: 1fr;
        }
        .content {
            padding: 20px;
        }
    }

    @media (max-width: 576px) {
        .sidebar {
            width: 60px;
            padding: 10px 0;
        }
        .sidebar-header {
            padding: 10px;
        }
        .sidebar-title,
        .sidebar-version {
            display: none;
        }
        .sidebar-nav-link {
            padding: 12px 10px;
            font-size: 0.8em;
            text-align: center;
        }
        .main-wrapper {
            margin-left: 60px;
        }
        .header {
            padding: 20px;
        }
        .content {
            padding: 15px;
        }
    }

    /* Subtab Navigation */
    .subtabs {
        display: flex;
        flex-wrap: wrap;
        gap: 5px;
        margin-bottom: 20px;
        border-bottom: 2px solid var(--border-color);
        padding-bottom: 0;
    }

    .subtab-button {
        padding: 10px 16px;
        background-color: #f5f5f5;
        border: 1px solid var(--border-color);
        border-bottom: none;
        border-radius: 6px 6px 0 0;
        cursor: pointer;
        font-size: 0.9em;
        color: #666;
        transition: all 0.2s ease;
    }

    .subtab-button:hover {
        background-color: #e8e8e8;
        color: var(--primary-color);
    }

    .subtab-button.active {
        background-color: white;
        color: var(--primary-color);
        font-weight: 600;
        border-bottom: 2px solid white;
        margin-bottom: -2px;
    }

    .subtab-content {
        display: none;
        animation: fadeIn 0.3s ease-in;
    }

    .subtab-content.active {
        display: block;
    }

    @keyframes fadeIn {
        from { opacity: 0; }
        to { opacity: 1; }
    }

    /* Visualization Styles */
    .viz-container {
        margin: 30px 0;
        padding: 20px;
        background-color: #f9f9f9;
        border-radius: 8px;
        border: 1px solid var(--border-color);
    }

    .viz-container h3 {
        margin-top: 0;
        color: var(--primary-color);
        font-size: 1.3em;
        margin-bottom: 10px;
    }

    .viz-description {
        color: #666;
        font-size: 0.95em;
        margin-bottom: 15px;
        font-style: italic;
    }

    .viz-image {
        width: 65%;
        max-width: 80%;
        height: auto;
        display: block;
        margin: 0 auto;
        border: 1px solid #ddd;
        border-radius: 4px;
        background-color: white;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        cursor: zoom-in;
    }

    .viz-image:hover {
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }

    /* Interactive Plot Controls */
    .plot-mode-toggle {
        display: flex;
        gap: 10px;
        margin-bottom: 20px;
        justify-content: center;
    }

    .toggle-btn {
        padding: 10px 24px;
        border: 2px solid var(--primary-color);
        background: white;
        color: var(--primary-color);
        font-weight: 600;
        border-radius: 6px;
        cursor: pointer;
        transition: all 0.3s ease;
        font-size: 14px;
    }

    .toggle-btn:hover {
        background: #f0f8ff;
    }

    .toggle-btn.active {
        background: var(--primary-color);
        color: white;
    }

    .plot-view {
        display: none;
    }

    .plot-view.active {
        display: block;
    }

    .plot-controls {
        background: #f8f9fa;
        border: 1px solid #dee2e6;
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
    }

    .controls-section h4 {
        margin-top: 0;
        margin-bottom: 16px;
        color: var(--text-color);
        font-size: 16px;
        font-weight: 600;
    }

    .control-group {
        margin-bottom: 20px;
    }

    .control-group label {
        display: block;
        font-weight: 600;
        margin-bottom: 8px;
        color: var(--text-color);
    }

    .controls-flex-container {
        display: flex;
        gap: 20px;
        margin-bottom: 20px;
    }

    .controls-flex-container .control-group {
        flex: 1;
        margin-bottom: 0;
    }

    .threshold-slider {
        width: 100%;
        height: 6px;
        border-radius: 3px;
        background: #ddd;
        outline: none;
        -webkit-appearance: none;
    }

    .threshold-slider::-webkit-slider-thumb {
        -webkit-appearance: none;
        appearance: none;
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: var(--primary-color);
        cursor: pointer;
        transition: all 0.2s ease;
    }

    .threshold-slider::-webkit-slider-thumb:hover {
        transform: scale(1.2);
        background: #0056b3;
    }

    .threshold-slider::-moz-range-thumb {
        width: 18px;
        height: 18px;
        border-radius: 50%;
        background: var(--primary-color);
        cursor: pointer;
        border: none;
    }

    .checkbox-controls {
        display: flex;
        gap: 8px;
        margin-bottom: 12px;
    }

    .genotype-checkboxes {
        max-height: 250px;
        overflow-y: auto;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        padding: 12px;
        background: white;
    }

    .genotype-checkbox-item {
        display: flex;
        align-items: center;
        padding: 6px 0;
        border-bottom: 1px solid #f0f0f0;
    }

    .genotype-checkbox-item:last-child {
        border-bottom: none;
    }

    .genotype-checkbox-item input[type="checkbox"] {
        margin-right: 10px;
        width: 18px;
        height: 18px;
        cursor: pointer;
    }

    .genotype-checkbox-item label {
        margin: 0;
        font-weight: normal;
        cursor: pointer;
        flex: 1;
        display: flex;
        justify-content: space-between;
        align-items: center;
    }

    .genotype-sample-count {
        color: #666;
        font-size: 0.9em;
        margin-left: 8px;
    }

    .btn-small {
        padding: 6px 12px;
        font-size: 13px;
        border: 1px solid #dee2e6;
        background: white;
        color: var(--text-color);
        border-radius: 4px;
        cursor: pointer;
        transition: all 0.2s ease;
    }

    .btn-small:hover {
        background: #e9ecef;
        border-color: #adb5bd;
    }

    .btn-primary {
        padding: 10px 20px;
        background: var(--primary-color);
        color: white;
        border: none;
        border-radius: 6px;
        font-weight: 600;
        cursor: pointer;
        transition: all 0.2s ease;
        font-size: 14px;
        width: 100%;
    }

    .btn-primary:hover {
        background: #0056b3;
        transform: translateY(-1px);
        box-shadow: 0 4px 8px rgba(0,0,0,0.15);
    }

    .plotly-chart {
        min-height: 500px;
        border: 1px solid #ddd;
        border-radius: 4px;
        background: white;
    }

    /* Summary Statistics Panel */
    .stats-panel {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        border-radius: 8px;
        padding: 20px;
        margin-bottom: 20px;
        box-shadow: 0 4px 6px rgba(0,0,0,0.1);
    }

    .stats-grid {
        display: grid;
        grid-template-columns: repeat(3, 1fr);
        gap: 20px;
    }

    .stat-item {
        text-align: center;
        background: rgba(255, 255, 255, 0.15);
        padding: 15px;
        border-radius: 6px;
        backdrop-filter: blur(10px);
    }

    .stat-label {
        font-size: 12px;
        font-weight: 600;
        text-transform: uppercase;
        letter-spacing: 0.5px;
        color: rgba(255, 255, 255, 0.9);
        margin-bottom: 8px;
    }

    .stat-value {
        font-size: 28px;
        font-weight: 700;
        color: white;
        line-height: 1;
    }

    @media (max-width: 768px) {
        .stats-grid {
            grid-template-columns: 1fr;
            gap: 12px;
        }
    }

    /* Download Options */
    .download-options {
        display: flex;
        gap: 8px;
        margin-bottom: 12px;
        flex-wrap: wrap;
    }

    .btn-download {
        flex: 1;
        min-width: 100px;
        padding: 10px 16px;
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        border: none;
        border-radius: 6px;
        font-weight: 600;
        cursor: pointer;
        transition: all 0.3s ease;
        font-size: 13px;
    }

    .btn-download:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 12px rgba(102, 126, 234, 0.4);
    }

    .resolution-control {
        margin-top: 8px;
    }

    .small-label {
        font-size: 13px;
        font-weight: 600;
        display: block;
        margin-bottom: 6px;
        color: var(--text-color);
    }

    .resolution-select {
        width: 100%;
        padding: 8px 12px;
        border: 1px solid #dee2e6;
        border-radius: 4px;
        background: white;
        font-size: 13px;
        cursor: pointer;
        transition: border-color 0.2s ease;
    }

    .resolution-select:hover {
        border-color: var(--primary-color);
    }

    .resolution-select:focus {
        outline: none;
        border-color: var(--primary-color);
        box-shadow: 0 0 0 3px rgba(0, 123, 255, 0.1);
    }

    @media print {
        .viz-container {
            page-break-inside: avoid;
        }
        .subtab-button:not(.active) {
            display: none;
        }
        .subtab-content {
            display: block !important;
        }
        .plot-mode-toggle,
        .plot-controls {
            display: none !important;
        }
        .plot-view {
            display: none !important;
        }
        .plot-view.active {
            display: block !important;
        }
    }

    /* Species Composition Styling */
    .stats-summary {
        display: flex;
        gap: 20px;
        margin-bottom: 30px;
        flex-wrap: wrap;
    }

    .stat-box {
        flex: 1;
        min-width: 200px;
        padding: 20px;
        background: #f8f9fa;
        border-radius: 8px;
        text-align: center;
        border-left: 4px solid #007bff;
    }

    .stat-box.warning {
        border-left-color: #ff9800;
        background: #fff3e0;
    }

    .stat-box.alert {
        border-left-color: #f44336;
        background: #ffebee;
    }

    .stat-box strong {
        font-size: 2em;
        display: block;
        margin-bottom: 5px;
        color: var(--text-color);
    }

    .table-container {
        overflow-x: auto;
        margin-top: 20px;
    }

    .data-table {
        width: 100%;
        border-collapse: collapse;
        background: white;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    }

    .data-table th {
        background: #007bff;
        color: white;
        padding: 12px;
        text-align: left;
        cursor: pointer;
        font-weight: 600;
        position: sticky;
        top: 0;
        z-index: 10;
    }

    .data-table th:hover {
        background: #0056b3;
    }

    .data-table td {
        padding: 10px 12px;
        border-bottom: 1px solid #ddd;
    }

    .data-table tbody tr:hover {
        background: #e3f2fd;
    }

    .data-table tbody tr.alert-row {
        background-color: #ffebee;
    }

    .data-table tbody tr.alert-row:hover {
        background-color: #ffcdd2;
    }

    .data-table tbody tr.warning-row {
        background-color: #fff3e0;
    }

    .data-table tbody tr.warning-row:hover {
        background-color: #ffe0b2;
    }

    .badge {
        display: inline-block;
        padding: 4px 8px;
        border-radius: 4px;
        font-size: 0.85em;
        font-weight: 600;
        margin: 2px;
    }

    .badge-warning {
        background: #ff9800;
        color: white;
    }

    .badge-alert {
        background: #f44336;
        color: white;
    }

    .composition-cell {
        font-family: 'Courier New', monospace;
        font-size: 0.9em;
        max-width: 400px;
        word-wrap: break-word;
    }

    .filter-container {
        margin: 20px 0;
        padding: 15px;
        background: #f8f9fa;
        border-radius: 4px;
    }

    .filter-container label {
        font-weight: 600;
        margin-right: 10px;
    }

    .filter-container select {
        padding: 8px 12px;
        border: 1px solid #ccc;
        border-radius: 4px;
        font-size: 1em;
        cursor: pointer;
    }
</style>
"""

# HTML Report Template
HTML_REPORT_TEMPLATE = """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <meta name="generator" content="BOLDGenotyper">
    <title>{{ organism }} - BOLDGenotyper Report</title>
    {{ css | safe }}
    <!-- Plotly.js for interactive charts -->
    <script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
</head>
<body>
    <div class="container">
        <!-- Left Sidebar Navigation -->
        <aside class="sidebar">
            <div class="sidebar-header">
                <div class="sidebar-title">{{ organism }}</div>
                <div class="sidebar-version">BOLDGenotyper v{{ version }}</div>
            </div>
            <nav>
                <ul class="sidebar-nav">
                    <li class="sidebar-nav-item">
                        <a href="#section-executive-summary" class="sidebar-nav-link">Summary</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-parameters" class="sidebar-nav-link">Parameters</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-visualizations" class="sidebar-nav-link">Visualizations</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-methods" class="sidebar-nav-link">Methods</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-assignment" class="sidebar-nav-link">Assignment</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-quality-control" class="sidebar-nav-link">Quality Control</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-taxonomy" class="sidebar-nav-link">Taxonomy</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-species-composition" class="sidebar-nav-link">Species Composition</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-geographic" class="sidebar-nav-link">Geography</a>
                    </li>
                    <li class="sidebar-nav-item">
                        <a href="#section-metadata-summary" class="sidebar-nav-link">Metadata Summary</a>
                    </li>
                </ul>
            </nav>
        </aside>

        <!-- Main Content Area -->
        <div class="main-wrapper">
            <!-- Header -->
            <div class="header">
                <h1>{{ organism }}</h1>
                <div class="subtitle">BOLDGenotyper Analysis Report</div>
                <div class="timestamp">Generated: {{ timestamp }}</div>
            </div>

            <!-- Quick Stats -->
            {% if quick_stats %}
            <div class="quick-stats">
                {% for stat in quick_stats %}
                <div class="quick-stat">
                    <div class="number">{{ stat.value }}</div>
                    <div class="label">{{ stat.label }}</div>
                </div>
                {% endfor %}
            </div>
            {% endif %}

            <!-- Main Content -->
            <div class="content">
                {{ content | safe }}
            </div>

            <!-- Footer -->
            <div class="footer">
                <p>
                    <strong>BOLDGenotyper</strong> - Automated genotyping pipeline for BOLD barcode data<br>
                    <a href="https://github.com/SymbioSeas/BOLDGenotyper" target="_blank">GitHub Repository</a> | Created by Steph Smith
                </p>
                <p style="margin-top: 15px; font-size: 0.85em;">
                    Report generated on {{ timestamp }}<br>
                    For questions or issues, please visit the GitHub repository.
                </p>
            </div>
        </div>
    </div>

    <script>
        // Switch between subtabs
        function switchSubtab(tabId) {
            // Find the parent section to scope the tab switching
            const tabElement = document.getElementById(tabId);
            if (!tabElement) return;

            const parentSection = tabElement.closest('.section');
            if (!parentSection) return;

            // Hide all subtab contents in this section
            const allContents = parentSection.querySelectorAll('.subtab-content');
            allContents.forEach(content => {
                content.classList.remove('active');
            });

            // Deactivate all subtab buttons in this section
            const allButtons = parentSection.querySelectorAll('.subtab-button');
            allButtons.forEach(button => {
                button.classList.remove('active');
            });

            // Show selected tab content
            tabElement.classList.add('active');

            // Activate corresponding button
            const buttonIndex = Array.from(allContents).indexOf(tabElement);
            if (buttonIndex >= 0 && allButtons[buttonIndex]) {
                allButtons[buttonIndex].classList.add('active');
            }
        }

        // Highlight active section in sidebar based on scroll position
        document.addEventListener('DOMContentLoaded', function() {
            const sections = document.querySelectorAll('.section[id^="section-"]');
            const navLinks = document.querySelectorAll('.sidebar-nav-link');

            function updateActiveLink() {
                let currentSection = '';
                sections.forEach(section => {
                    const sectionTop = section.offsetTop - 100;
                    if (window.pageYOffset >= sectionTop) {
                        currentSection = section.getAttribute('id');
                    }
                });

                navLinks.forEach(link => {
                    link.classList.remove('active');
                    if (link.getAttribute('href') === '#' + currentSection) {
                        link.classList.add('active');
                    }
                });
            }

            window.addEventListener('scroll', updateActiveLink);
            updateActiveLink(); // Initial call

            // Smooth scrolling for anchor links
            navLinks.forEach(link => {
                link.addEventListener('click', function(e) {
                    e.preventDefault();
                    const targetId = this.getAttribute('href').substring(1);
                    const targetElement = document.getElementById(targetId);
                    if (targetElement) {
                        window.scrollTo({
                            top: targetElement.offsetTop - 20,
                            behavior: 'smooth'
                        });
                    }
                });
            });

            // Initialize interactive plots when DOM is loaded
            initializeInteractivePlots();
        });

        // Global state for interactive plots
        const plotState = {};

        // Global color state - shared across all plots
        const globalColors = {};
        let allGenotypes = new Set();

        // Initialize all interactive plots
        function initializeInteractivePlots() {
            // Find all embedded plot data
            const plotDataElements = document.querySelectorAll('script[type="application/json"][id^="plot-data-"]');

            plotDataElements.forEach(dataElement => {
                const idx = parseInt(dataElement.id.replace('plot-data-', ''));
                const plotData = JSON.parse(dataElement.textContent);

                // Collect all genotypes across all plots
                if (plotData.genotypes) {
                    plotData.genotypes.forEach(g => allGenotypes.add(g));
                }

                // Initialize global colors from first plot's colors if not already set
                if (plotData.colors && Object.keys(globalColors).length === 0) {
                    Object.assign(globalColors, plotData.colors);
                }

                // Initialize plot state
                plotState[idx] = {
                    data: plotData,
                    threshold: 0,
                    selectedGenotypes: new Set(plotData.genotypes || []),
                    selectedBasins: new Set(plotData.ocean_basins || [])
                };

                // Populate genotype checkboxes
                populateGenotypeCheckboxes(idx);

                // Populate ocean basin checkboxes if available
                if (plotData.ocean_basins && plotData.ocean_basins.length > 0) {
                    populateBasinCheckboxes(idx);
                }
            });

            // Initialize global color picker
            initializeGlobalColorPicker();
        }

        // Initialize global color picker panel
        function initializeGlobalColorPicker() {
            const container = document.getElementById('global-color-picker');
            if (!container || allGenotypes.size === 0) return;

            let html = '<h4>Global Genotype Colors</h4>';
            html += '<p style="font-size: 0.85em; color: #666; margin-bottom: 10px;">Adjust colors for all interactive plots simultaneously</p>';
            html += '<div style="max-height: 300px; overflow-y: auto;">';

            const sortedGenotypes = Array.from(allGenotypes).sort();
            sortedGenotypes.forEach(genotype => {
                const color = globalColors[genotype] || '#888888';
                // Escape genotype name for use in HTML IDs and JavaScript
                const escapedId = genotype.replace(/[^a-zA-Z0-9-_]/g, '_');
                html += '<div style="display: flex; align-items: center; margin-bottom: 8px; padding: 5px; border-radius: 4px; background: #f8f9fa;">';
                html += `<input type="color" id="color-${escapedId}" value="${color}" `;
                html += `onchange="updateGlobalColor(\\`${genotype}\\`, this.value)" `;
                html += 'style="width: 40px; height: 30px; border: 1px solid #ccc; border-radius: 4px; cursor: pointer; margin-right: 10px;">';
                html += `<label for="color-${escapedId}" style="flex: 1; font-size: 0.9em; cursor: pointer;">${genotype}</label>`;
                html += '</div>';
            });

            html += '</div>';
            html += '<button class="btn" onclick="resetAllColors()" style="margin-top: 10px; width: 100%;">Reset to Default Colors</button>';

            container.innerHTML = html;
        }

        // Update a single genotype's global color
        function updateGlobalColor(genotype, color) {
            globalColors[genotype] = color;

            // Update all plots
            Object.keys(plotState).forEach(idx => {
                updatePlot(parseInt(idx));
            });
        }

        // Reset all colors to defaults
        function resetAllColors() {
            // Reset globalColors to original values from plot data
            Object.keys(plotState).forEach(idx => {
                const plotData = plotState[idx].data;
                if (plotData.colors) {
                    Object.assign(globalColors, plotData.colors);
                }
            });

            // Update color picker inputs
            Object.keys(globalColors).forEach(genotype => {
                const escapedId = genotype.replace(/[^a-zA-Z0-9-_]/g, '_');
                const input = document.getElementById(`color-${escapedId}`);
                if (input) {
                    input.value = globalColors[genotype];
                }
            });

            // Re-render all plots
            Object.keys(plotState).forEach(idx => {
                updatePlot(parseInt(idx));
            });
        }

        // Toggle between static and interactive plot views
        function showStaticPlot(plotId, idx) {
            // Update button states
            const buttons = document.querySelectorAll(`button[onclick*="'${plotId}'"]`);
            buttons.forEach(btn => {
                if (btn.textContent.includes('Static')) {
                    btn.classList.add('active');
                } else {
                    btn.classList.remove('active');
                }
            });

            // Show static, hide interactive
            document.getElementById(`${plotId}-static`).style.display = 'block';
            document.getElementById(`${plotId}-interactive`).style.display = 'none';
        }

        function showInteractivePlot(plotId, idx) {
            // Update button states
            const buttons = document.querySelectorAll(`button[onclick*="'${plotId}'"]`);
            buttons.forEach(btn => {
                if (btn.textContent.includes('Interactive')) {
                    btn.classList.add('active');
                } else {
                    btn.classList.remove('active');
                }
            });

            // Show interactive, hide static
            document.getElementById(`${plotId}-static`).style.display = 'none';
            document.getElementById(`${plotId}-interactive`).style.display = 'block';

            // Render the plot if not already rendered
            if (!plotState[idx].rendered) {
                updatePlot(idx);
                plotState[idx].rendered = true;
            }
        }

        // Populate genotype checkboxes
        function populateGenotypeCheckboxes(idx) {
            const container = document.getElementById(`genotype-checkboxes-${idx}`);
            if (!container) return;

            const state = plotState[idx];
            const plotData = state.data;

            container.innerHTML = '';

            const genotypes = plotData.genotypes || [];
            const sampleCounts = plotData.sample_counts || {};

            genotypes.forEach(genotype => {
                const div = document.createElement('div');
                div.className = 'genotype-checkbox-item';

                const checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.id = `genotype-${idx}-${genotype}`;
                checkbox.value = genotype;
                checkbox.checked = state.selectedGenotypes.has(genotype);
                checkbox.onchange = () => {
                    toggleGenotype(idx, genotype);
                    updatePlot(idx);
                };

                const label = document.createElement('label');
                label.htmlFor = checkbox.id;
                label.innerHTML = `
                    <span>${genotype}</span>
                    <span class="genotype-sample-count">n=${sampleCounts[genotype] || 0}</span>
                `;

                div.appendChild(checkbox);
                div.appendChild(label);
                container.appendChild(div);
            });
        }

        // Toggle genotype selection
        function toggleGenotype(idx, genotype) {
            const state = plotState[idx];
            if (state.selectedGenotypes.has(genotype)) {
                state.selectedGenotypes.delete(genotype);
            } else {
                state.selectedGenotypes.add(genotype);
            }
        }

        // Select all genotypes
        function selectAllGenotypes(idx) {
            const state = plotState[idx];
            const genotypes = state.data.genotypes || [];

            genotypes.forEach(g => state.selectedGenotypes.add(g));

            // Update checkbox states
            genotypes.forEach(g => {
                const checkbox = document.getElementById(`genotype-${idx}-${g}`);
                if (checkbox) checkbox.checked = true;
            });

            updatePlot(idx);
        }

        // Deselect all genotypes
        function deselectAllGenotypes(idx) {
            const state = plotState[idx];

            state.selectedGenotypes.clear();

            // Update checkbox states
            const checkboxes = document.querySelectorAll(`input[id^="genotype-${idx}-"]`);
            checkboxes.forEach(cb => cb.checked = false);

            updatePlot(idx);
        }

        // Populate ocean basin checkboxes
        function populateBasinCheckboxes(idx) {
            const container = document.getElementById(`basin-checkboxes-${idx}`);
            if (!container) return;

            const state = plotState[idx];
            const plotData = state.data;
            const basins = plotData.ocean_basins || [];

            container.innerHTML = '';

            basins.forEach(basin => {
                const div = document.createElement('div');
                div.className = 'genotype-checkbox-item';

                const checkbox = document.createElement('input');
                checkbox.type = 'checkbox';
                checkbox.id = `basin-${idx}-${basin}`;
                checkbox.value = basin;
                checkbox.checked = state.selectedBasins.has(basin);
                checkbox.onchange = () => {
                    toggleBasin(idx, basin);
                    updatePlot(idx);
                };

                const label = document.createElement('label');
                label.htmlFor = checkbox.id;
                label.textContent = basin;

                div.appendChild(checkbox);
                div.appendChild(label);
                container.appendChild(div);
            });
        }

        // Toggle ocean basin selection
        function toggleBasin(idx, basin) {
            const state = plotState[idx];
            if (state.selectedBasins.has(basin)) {
                state.selectedBasins.delete(basin);
            } else {
                state.selectedBasins.add(basin);
            }
        }

        // Select all ocean basins
        function selectAllBasins(idx) {
            const state = plotState[idx];
            const basins = state.data.ocean_basins || [];

            basins.forEach(b => state.selectedBasins.add(b));

            // Update checkbox states
            basins.forEach(b => {
                const checkbox = document.getElementById(`basin-${idx}-${b}`);
                if (checkbox) checkbox.checked = true;
            });

            updatePlot(idx);
        }

        // Deselect all ocean basins
        function deselectAllBasins(idx) {
            const state = plotState[idx];

            state.selectedBasins.clear();

            // Update checkbox states
            const checkboxes = document.querySelectorAll(`input[id^="basin-${idx}-"]`);
            checkboxes.forEach(cb => cb.checked = false);

            updatePlot(idx);
        }

        // Update threshold value
        function updateThreshold(idx) {
            const slider = document.getElementById(`threshold-${idx}`);
            const valueDisplay = document.getElementById(`threshold-value-${idx}`);

            if (slider && valueDisplay) {
                plotState[idx].threshold = parseInt(slider.value);
                valueDisplay.textContent = slider.value;
            }
        }

        // Update summary statistics
        function updateStatistics(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // Calculate filtered genotypes
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            // Calculate total samples in filtered data
            let totalSamples = 0;
            filteredGenotypes.forEach(g => {
                totalSamples += plotData.sample_counts[g] || 0;
            });

            // Calculate percentage of total data
            let allSamples = 0;
            (plotData.genotypes || []).forEach(g => {
                allSamples += plotData.sample_counts[g] || 0;
            });
            const percentage = allSamples > 0 ? (totalSamples / allSamples * 100).toFixed(1) : 0;

            // Update DOM elements
            const genotypesStat = document.getElementById(`stat-genotypes-${idx}`);
            const samplesStat = document.getElementById(`stat-samples-${idx}`);
            const percentageStat = document.getElementById(`stat-percentage-${idx}`);

            if (genotypesStat) {
                genotypesStat.textContent = `${filteredGenotypes.length} / ${plotData.genotypes.length}`;
            }
            if (samplesStat) {
                samplesStat.textContent = totalSamples.toLocaleString();
            }
            if (percentageStat) {
                percentageStat.textContent = `${percentage}%`;
            }
        }

        // Update/render the Plotly chart
        function updatePlot(idx) {
            const state = plotState[idx];
            const plotData = state.data;
            const plotType = plotData.plot_type;

            // Update statistics
            updateStatistics(idx);

            if (plotType === 'faceted_abundance') {
                renderFacetedPlot(idx);
            } else if (plotType === 'distribution_map') {
                renderDistributionMap(idx);
            } else {
                renderStackedBarPlot(idx);
            }
        }

        // Render stacked bar plot (for relative and total abundance)
        function renderStackedBarPlot(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // Filter genotypes based on selection and threshold
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            if (filteredGenotypes.length === 0) {
                // Show message if no genotypes selected
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No genotypes selected or meet the minimum sample threshold.</div>';
                }
                return;
            }

            // Prepare Plotly traces (one per genotype)
            const traces = [];
            const basins = plotData.basins || [];
            const basinLabels = plotData.basin_labels || {};
            const colors = plotData.colors || [];
            const counts = plotData.counts || {};

            filteredGenotypes.forEach((genotype, i) => {
                const genotypeCounts = counts[genotype] || {};
                const y = basins.map(basin => genotypeCounts[basin] || 0);

                traces.push({
                    name: genotype,
                    x: basins.map(b => basinLabels[b] || b),
                    y: y,
                    type: 'bar',
                    marker: {
                        color: globalColors[genotype] || colors[genotype] || '#888',
                        line: {
                            color: 'black',
                            width: 1
                        }
                    }
                });
            });

            // Layout configuration
            const isRelative = plotData.plot_type === 'relative_abundance';
            const layout = {
                barmode: 'stack',
                title: '',
                xaxis: {
                    title: 'Ocean Basin',
                    tickangle: 0
                },
                yaxis: {
                    title: isRelative ? 'Relative Abundance' : 'Total Sample Count',
                    rangemode: 'tozero'
                },
                legend: {
                    orientation: 'v',
                    x: 1.02,
                    y: 1,
                    xanchor: 'left'
                },
                margin: {
                    l: 60,
                    r: 200,
                    t: 40,
                    b: 80
                },
                height: 300,
                hovermode: 'closest'
            };

            if (isRelative) {
                layout.barnorm = 'fraction';
                layout.yaxis.tickformat = '.0%';
            }

            // Render plot
            const config = {
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            };

            Plotly.newPlot(`plot-${idx}`, traces, layout, config);
        }

        // Render faceted plot
        function renderFacetedPlot(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // For faceted plots, we show all facets but filter genotypes
            const facets = plotData.facets || {};
            const facetNames = Object.keys(facets);

            if (facetNames.length === 0) {
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No facets available.</div>';
                }
                return;
            }

            // Create traces for all facets
            const traces = [];
            const basins = plotData.basins || [];
            const basinLabels = plotData.basin_labels || {};
            const colors = plotData.colors || [];

            facetNames.forEach(facetName => {
                const facet = facets[facetName];
                const genotypes = facet.genotypes || [];

                genotypes.forEach(genotype => {
                    // Check if this genotype is selected
                    if (!state.selectedGenotypes.has(facetName)) return;

                    const genotypeCounts = facet.counts[genotype] || {};
                    const y = basins.map(basin => genotypeCounts[basin] || 0);

                    traces.push({
                        name: facetName,
                        x: basins.map(b => basinLabels[b] || b),
                        y: y,
                        type: 'bar',
                        marker: {
                            color: globalColors[genotype] || colors[genotype] || '#888',
                            line: {
                                color: 'black',
                                width: 1
                            }
                        }
                    });
                });
            });

            if (traces.length === 0) {
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No genotypes selected.</div>';
                }
                return;
            }

            // Layout configuration
            const layout = {
                barmode: 'stack',
                title: '',
                xaxis: {
                    title: 'Ocean Basin',
                    tickangle: 0
                },
                yaxis: {
                    title: 'Sample Count',
                    rangemode: 'tozero'
                },
                legend: {
                    orientation: 'v',
                    x: 1.02,
                    y: 1,
                    xanchor: 'left'
                },
                margin: {
                    l: 60,
                    r: 200,
                    t: 40,
                    b: 80
                },
                height: 450,
                hovermode: 'closest'
            };

            // Render plot
            const config = {
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            };

            Plotly.newPlot(`plot-${idx}`, traces, layout, config);
        }

        // Render distribution map
        function renderDistributionMap(idx) {
            const state = plotState[idx];
            const plotData = state.data;

            // Filter genotypes based on selection and threshold
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            if (filteredGenotypes.length === 0) {
                const plotDiv = document.getElementById(`plot-${idx}`);
                if (plotDiv) {
                    plotDiv.innerHTML = '<div style="padding: 40px; text-align: center; color: #666;">No genotypes selected or meet the minimum sample threshold.</div>';
                }
                return;
            }

            // Prepare Plotly traces (one per genotype)
            const traces = [];
            const locations = plotData.locations || {};
            const colors = plotData.colors || {};

            filteredGenotypes.forEach(genotype => {
                const genoLocations = locations[genotype];
                if (!genoLocations) return;

                let lats = genoLocations.latitudes || [];
                let lons = genoLocations.longitudes || [];
                let sizes = genoLocations.sizes || [];
                const oceanBasins = genoLocations.ocean_basins || [];

                // Filter by ocean basin if basin data is available and basins are selected
                if (oceanBasins.length > 0 && state.selectedBasins.size > 0) {
                    const filteredIndices = [];
                    oceanBasins.forEach((basin, i) => {
                        if (state.selectedBasins.has(basin)) {
                            filteredIndices.push(i);
                        }
                    });

                    // Filter all arrays by selected basins
                    lats = filteredIndices.map(i => lats[i]);
                    lons = filteredIndices.map(i => lons[i]);
                    sizes = filteredIndices.map(i => sizes[i]);
                }

                // Skip if no points remaining after filtering
                if (lats.length === 0) return;

                // Scale marker sizes for better visibility
                const scaledSizes = sizes.map(s => Math.max(5, Math.min(20, s * 3)));

                traces.push({
                    type: 'scattergeo',
                    mode: 'markers',
                    name: genotype,
                    lat: lats,
                    lon: lons,
                    marker: {
                        size: scaledSizes,
                        color: globalColors[genotype] || colors[genotype] || '#888',
                        line: {
                            color: 'black',
                            width: 0.5
                        },
                        opacity: 0.85
                    },
                    hovertemplate: '<b>%{text}</b><br>Lat: %{lat:.2f}<br>Lon: %{lon:.2f}<extra></extra>',
                    text: lats.map(() => genotype)
                });
            });

            // Layout configuration
            const layout = {
                title: '',
                geo: {
                    projection: {
                        type: 'natural earth'
                    },
                    showland: true,
                    landcolor: 'rgb(243, 243, 243)',
                    coastlinecolor: 'rgb(204, 204, 204)',
                    showocean: true,
                    oceancolor: 'rgb(230, 245, 255)',
                    showcountries: true,
                    countrycolor: 'rgb(204, 204, 204)',
                    showlakes: true,
                    lakecolor: 'rgb(230, 245, 255)',
                    resolution: 50
                },
                legend: {
                    orientation: 'v',
                    x: 1.02,
                    y: 1,
                    xanchor: 'left'
                },
                margin: {
                    l: 0,
                    r: 200,
                    t: 40,
                    b: 0
                },
                height: 450,
                hovermode: 'closest'
            };

            // Render plot
            const config = {
                responsive: true,
                displayModeBar: true,
                displaylogo: false,
                modeBarButtonsToRemove: ['lasso2d', 'select2d']
            };

            Plotly.newPlot(`plot-${idx}`, traces, layout, config);
        }

        // Download plot as image (PNG or SVG)
        function downloadPlot(idx, format) {
            const plotDiv = document.getElementById(`plot-${idx}`);
            if (!plotDiv) return;

            // Get selected resolution
            const resolutionSelect = document.getElementById(`resolution-${idx}`);
            const resolution = resolutionSelect ? resolutionSelect.value : '1200x900';
            const [width, height] = resolution.split('x').map(Number);

            // Determine filename based on plot type
            const state = plotState[idx];
            const plotType = state.data.plot_type || 'plot';
            const filename = `${plotType}_filtered_${idx}`;

            // Use Plotly's built-in download functionality
            Plotly.downloadImage(plotDiv, {
                format: format,
                width: width,
                height: height,
                filename: filename
            });
        }

        // Download filtered data as CSV
        function downloadCSV(idx) {
            const state = plotState[idx];
            const plotData = state.data;
            const plotType = plotData.plot_type;

            // Filter genotypes based on current selection and threshold
            const filteredGenotypes = (plotData.genotypes || []).filter(g => {
                if (!state.selectedGenotypes.has(g)) return false;
                const sampleCount = plotData.sample_counts[g] || 0;
                return sampleCount >= state.threshold;
            });

            if (filteredGenotypes.length === 0) {
                alert('No data to export with current filters');
                return;
            }

            let csvContent = '';

            if (plotType === 'distribution_map') {
                // CSV for distribution map: genotype, latitude, longitude, sample_count
                csvContent = 'Genotype,Latitude,Longitude,Sample_Count\\n';

                filteredGenotypes.forEach(genotype => {
                    const locations = plotData.locations[genotype];
                    if (!locations) return;

                    const lats = locations.latitudes || [];
                    const lons = locations.longitudes || [];
                    const sizes = locations.sizes || [];

                    for (let i = 0; i < lats.length; i++) {
                        csvContent += `"${genotype}",${lats[i]},${lons[i]},${sizes[i]}\n`;
                    }
                });
            } else {
                // CSV for bar plots: genotype, basin, count, sample_total
                csvContent = 'Genotype,Ocean_Basin,Count,Total_Samples\\n';

                const basins = plotData.basins || [];
                const counts = plotData.counts || {};

                filteredGenotypes.forEach(genotype => {
                    const genotypeCounts = counts[genotype] || {};
                    const totalSamples = plotData.sample_counts[genotype] || 0;

                    basins.forEach(basin => {
                        const count = genotypeCounts[basin] || 0;
                        if (count > 0) {
                            csvContent += `"${genotype}","${basin}",${count},${totalSamples}\n`;
                        }
                    });
                });
            }

            // Create and download CSV file
            const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
            const link = document.createElement('a');
            const url = URL.createObjectURL(blob);

            link.setAttribute('href', url);
            link.setAttribute('download', `${plotType}_filtered_data_${idx}.csv`);
            link.style.visibility = 'hidden';
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
        }

        // Species composition table filtering
        function filterCompositionTable() {
            const filter = document.getElementById('comp-filter').value;
            const table = document.getElementById('composition-table');
            if (!table) return;

            const tbody = table.getElementsByTagName('tbody')[0];
            if (!tbody) return;

            const rows = tbody.getElementsByTagName('tr');

            for (let row of rows) {
                const isMulti = row.getAttribute('data-multi') === 'True';
                const isAmbiguous = row.getAttribute('data-ambiguous') === 'True';

                let show = false;
                if (filter === 'all') {
                    show = true;
                } else if (filter === 'multi' && isMulti) {
                    show = true;
                } else if (filter === 'ambiguous' && isAmbiguous) {
                    show = true;
                }

                row.style.display = show ? '' : 'none';
            }
        }
    </script>
</body>
</html>
"""


class HTMLReportBuilder:
    """
    Builder class for generating comprehensive HTML reports.

    Assembles pipeline results into a formatted HTML document with
    embedded tables, metrics, and visualizations.
    """

    def __init__(self, organism: str, version: Optional[str] = None):
        """
        Initialize HTML report builder.

        Parameters
        ----------
        organism : str
            Organism name
        version : str, optional
            BOLDGenotyper version (default: the installed package version)
        """
        if version is None:
            from . import __version__ as version
        self.organism = organism
        self.version = version
        self.sections = []
        self.quick_stats = []

    def add_quick_stat(self, value: str, label: str):
        """Add a quick stat to the header bar."""
        self.quick_stats.append({'value': value, 'label': label})

    def add_section(self, content: str):
        """Add a section to the report."""
        self.sections.append(content)

    def render(self) -> str:
        """
        Render the complete HTML report.

        Returns
        -------
        str
            Complete HTML document
        """
        try:
            from jinja2 import Template
        except ImportError:
            raise ImportError(
                "jinja2 is required for HTML report generation. "
                "Install it with: pip install jinja2"
            )

        from datetime import datetime

        template = Template(HTML_REPORT_TEMPLATE)

        return template.render(
            organism=self.organism,
            version=self.version,
            timestamp=datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            css=HTML_REPORT_CSS,
            quick_stats=self.quick_stats,
            content='\n'.join(self.sections)
        )


def _format_number(value: float, decimals: int = 2) -> str:
    """Format a number for display."""
    if pd.isna(value):
        return "N/A"
    if isinstance(value, (int, np.integer)):
        return f"{value:,}"
    return f"{value:,.{decimals}f}"


def _format_percentage(value: float, decimals: int = 1) -> str:
    """Format a percentage for display."""
    if pd.isna(value):
        return "N/A"
    return f"{value:.{decimals}f}%"


def _dataframe_to_html(df: pd.DataFrame, max_rows: int = 100) -> str:
    """
    Convert DataFrame to HTML table.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to convert
    max_rows : int, optional
        Maximum rows to display (default: 100)

    Returns
    -------
    str
        HTML table string
    """
    if df.empty:
        return '<p class="alert alert-info">No data available</p>'

    # Truncate if too long (count the rows first: the notice gives the true total)
    n_total = len(df)
    truncated = n_total > max_rows
    if truncated:
        df = df.head(max_rows)

    # Convert to HTML
    html = '<div class="table-container">\n'
    html += df.to_html(index=False, classes='', border=0, escape=True)
    html += '</div>\n'

    if truncated:
        html += (f'<p style="color: #888; font-size: 0.9em; margin-top: 10px;">'
                 f'Showing first {max_rows} rows of {n_total} total</p>\n')

    return html


def _compute_species_composition_from_df(
    df: pd.DataFrame,
    species_col: str = 'species',
    min_abundance_pct: float = 1.0
) -> pd.DataFrame:
    """Compute species composition per haplotype directly from annotated metadata."""
    required_cols = {'haplotype_id', species_col}
    if not required_cols.issubset(df.columns):
        return pd.DataFrame()

    assigned = df[df['haplotype_id'].notna()].copy()
    if assigned.empty:
        return pd.DataFrame()

    records = []
    for haplotype, group in assigned.groupby('haplotype_id'):
        species_series = group[species_col].fillna('Unknown')
        species_counts = species_series.value_counts()
        total = len(group)
        species_pcts = species_counts / total * 100
        significant = species_pcts[species_pcts >= min_abundance_pct]

        if significant.empty:
            continue

        comp_parts = [
            f"{sp}: {species_counts[sp]} ({species_pcts[sp]:.1f}%)"
            for sp in significant.index
        ]
        primary_species = species_counts.index[0]
        primary_pct = species_pcts.iloc[0]

        records.append({
            'haplotype_id': haplotype,
            'n_total_samples': total,
            'n_species': len(significant),
            'primary_species': primary_species,
            'primary_species_pct': f"{primary_pct:.1f}%",
            'species_composition': "; ".join(comp_parts),
            'is_multi_species': len(significant) > 1,
            'is_ambiguous': primary_pct < 70.0,
            'species_list': ", ".join(significant.index)
        })

    if not records:
        return pd.DataFrame()

    return pd.DataFrame(records).sort_values('n_total_samples', ascending=False)


def _build_executive_summary_section(
    assignment_summary: pd.DataFrame,
    annotated_df: pd.DataFrame
) -> str:
    """Build executive summary section with metric cards."""
    if assignment_summary.empty:
        return (
            '<div class="section" id="section-executive-summary">'
            '<h2>Analysis Summary</h2>'
            '<p class="alert alert-warning">No summary data available</p>'
            '</div>'
        )
    
    summary = assignment_summary.iloc[0] if len(assignment_summary) > 0 else {}

    html = '<div class="section" id="section-executive-summary">\n'
    html += '<h2>Analysis Summary</h2>\n'
    html += '<div class="metric-grid">\n'

    # Sample metrics
    html += '<div class="metric-card">\n'
    html += '  <div class="title">Unique Records</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_unique", 0), 0)}</div>\n'
    html += f'  <div class="description">{_format_number(summary.get("n_input", 0), 0)} input rows before duplicate removal</div>\n'
    html += '</div>\n'

    html += '<div class="metric-card">\n'
    html += '  <div class="title">QC-Passed</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_qc_pass", 0), 0)}</div>\n'
    pct = (summary.get("n_qc_pass", 0) / summary.get("n_unique", 1) * 100) if summary.get("n_unique", 0) > 0 else 0
    html += f'  <div class="description">{_format_percentage(pct)} of unique records</div>\n'
    html += '</div>\n'

    # Haplotype metrics
    html += '<div class="metric-card">\n'
    html += '  <div class="title">Haplotypes</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_haplotypes", 0), 0)}</div>\n'
    html += '  <div class="description">Unique haplotypes identified</div>\n'
    html += '</div>\n'

    html += '<div class="metric-card">\n'
    html += '  <div class="title">Species Assigned to Haplotypes</div>\n'
    html += f'  <div class="value">{_format_number(summary.get("n_species_haplotype_assigned", 0), 0)}</div>\n'
    html += '  <div class="description">Distinct species-level labels assigned to haplotypes</div>\n'
    html += '</div>\n'

    # Assignment metrics
    html += '<div class="metric-card">\n'
    html += '  <div class="title">Assigned to a Haplotype</div>\n'
    n_assigned_total = summary.get("n_assigned", 0) + summary.get("n_assigned_corrected", 0)
    html += f'  <div class="value">{_format_number(n_assigned_total, 0)}</div>\n'
    html += f'  <div class="description">{_format_percentage(summary.get("pct_assigned", 0))} of QC-passed records</div>\n'
    html += '</div>\n'

    # Geographic coverage
    if 'ocean_basin' in annotated_df.columns:
        n_basins = annotated_df['ocean_basin'].nunique()
        html += '<div class="metric-card">\n'
        html += '  <div class="title">Ocean Basins</div>\n'
        html += f'  <div class="value">{n_basins}</div>\n'
        html += '  <div class="description">Geographic coverage</div>\n'
        html += '</div>\n'

    html += '</div>\n'  # Close metric-grid
    html += '</div>\n'  # Close section

    return html


def _build_assignment_section(
    assignment_summary: pd.DataFrame,
    diagnostics_df: pd.DataFrame
) -> str:
    """Build genotype assignment results section."""
    html = '<div class="section" id="section-assignment">\n'
    html += '<h2>Genotype Assignment Results</h2>\n'
    html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 15px;">Full assignment data: <code>summary.csv</code> and <code>records.csv</code> (one row per record, with its status)</p>\n'

    if assignment_summary.empty:
        html += '<p class="alert alert-warning">No assignment summary data available</p>\n'
        html += '</div>\n'
        return html

    summary = assignment_summary.iloc[0]

    # Assignment Status Breakdown
    html += '<h3>Assignment Status Breakdown</h3>\n'
    html += '<div class="table-container">\n'
    html += '<table>\n'
    html += '<thead><tr><th>Status</th><th>Count</th><th>% of unique records</th></tr></thead>\n'
    html += '<tbody>\n'

    statuses = [
        ('n_assigned', 'assigned', 'badge-success'),
        ('n_assigned_corrected', 'assigned_corrected', 'badge-success'),
        ('n_excluded_error', 'excluded_error_singleton', 'badge-warning'),
        ('n_excluded_core', 'excluded_core_coverage', 'badge-warning'),
        ('n_failed_qc', 'failed_qc', 'badge-info'),
        ('n_no_sequence', 'no_sequence', 'badge-info'),
    ]
    n_unique = summary.get('n_unique', 0)

    for count_key, status_label, badge_class in statuses:
        count = summary.get(count_key, 0)
        pct = (count / n_unique * 100) if n_unique else 0
        html += f'<tr><td><span class="badge {badge_class}">{status_label}</span></td>'
        html += f'<td>{_format_number(count, 0)}</td>'
        html += f'<td>{_format_percentage(pct)}</td></tr>\n'

    html += '</tbody>\n</table>\n</div>\n'

    html += '</div>\n'  # Close section
    return html


def _build_taxonomy_section(
    taxonomy_dir: Path,
    organism: str
) -> str:
    """Build taxonomy section with consensus taxonomy and species composition."""
    html = '<div class="section" id="section-taxonomy">\n'
    html += '<h2>Taxonomy</h2>\n'

    # Load haplotype taxonomy
    consensus_tax_file = taxonomy_dir / f"{organism}_haplotype_taxonomy.csv"
    # Species composition is saved in haplotype_assignments/, not taxonomy/
    assignments_dir = taxonomy_dir.parent / 'haplotype_assignments'
    species_by_cons_file = assignments_dir / f"{organism}_species_composition.csv"
    # Fall back to taxonomy dir for backward compatibility
    if not species_by_cons_file.exists():
        species_by_cons_file = taxonomy_dir / f"{organism}_species_composition.csv"

    if consensus_tax_file.exists():
        html += '<h3>Haplotype Taxonomy</h3>\n'
        html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Full table: <code>haplotypes.csv</code></p>\n'
        try:
            tax_df = pd.read_csv(consensus_tax_file)
            # Select key columns
            display_cols = [c for c in ['haplotype_id', 'assigned_sp', 'assignment_level',
                                        'assignment_notes', 'majority_fraction'] if c in tax_df.columns]
            if display_cols:
                html += _dataframe_to_html(tax_df[display_cols], max_rows=50)
            else:
                html += _dataframe_to_html(tax_df, max_rows=50)
        except Exception as e:
            html += f'<p class="alert alert-warning">Could not load consensus taxonomy: {e}</p>\n'
    else:
        html += '<p class="alert alert-info">Consensus taxonomy file not found</p>\n'

    if species_by_cons_file.exists():
        html += '<h3>Species Composition per Haplotype</h3>\n'
        html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Full table: <code>haplotypes.csv</code> (<code>labels</code>, <code>majority_label</code>, <code>majority_fraction</code>)</p>\n'
        try:
            species_df = pd.read_csv(species_by_cons_file)
            html += _dataframe_to_html(species_df, max_rows=50)
        except Exception as e:
            html += f'<p class="alert alert-warning">Could not load species composition: {e}</p>\n'
    else:
        html += '<p class="alert alert-info">Species composition file not found</p>\n'

    html += '</div>\n'  # Close section
    return html


def _build_species_composition_section(
    results_dir: Path,
    organism: str,
    annotated_df: Optional[pd.DataFrame] = None
) -> str:
    """Build HTML section showing species composition of consensus groups."""
    html = '<div class="section" id="section-species-composition">\n'
    html += '<h2>Species Composition of Consensus Groups</h2>\n'

    # Load species composition file, or compute it if missing and data available
    comp_file = results_dir / f"{organism}_species_composition.csv"
    comp_df = pd.DataFrame()

    if comp_file.exists():
        try:
            comp_df = pd.read_csv(comp_file)
        except Exception as e:
            html += f'<p class="alert alert-warning">Could not load species composition file: {e}</p>\n'
            html += '</div>\n'
            return html
    elif annotated_df is not None and not annotated_df.empty:
        # Attempt to compute composition on the fly
        species_col = 'species' if 'species' in annotated_df.columns else 'species_name' if 'species_name' in annotated_df.columns else None
        if species_col:
            comp_df = _compute_species_composition_from_df(annotated_df, species_col=species_col)
            if not comp_df.empty:
                comp_file.parent.mkdir(parents=True, exist_ok=True)
                comp_df.to_csv(comp_file, index=False)
                html += '<p class="alert alert-info">Computed species composition from annotated data (file was missing in outputs).</p>\n'

    if comp_df.empty:
        html += '<p class="alert alert-info">Species composition analysis not available. This may be because:</p>\n'
        html += '<ul>\n'
        html += '<li>No species information was present in the input metadata</li>\n'
        html += '<li>The pipeline was run before this feature was added</li>\n'
        html += '</ul>\n'
        html += '</div>\n'
        return html

    # Add explanation
    html += '<p style="margin-bottom: 20px;">'
    html += 'This table shows which species names were reported for sequences within each haplotype. '
    html += 'When a haplotype\u2019s members carry more than one species label, the counts are shown; they do not decide which label is correct.'
    html += '</p>\n'

    # Summary statistics
    total_groups = len(comp_df)
    multi_species = comp_df['is_multi_species'].sum()
    ambiguous = comp_df['is_ambiguous'].sum()

    html += '<div class="stats-summary">\n'
    html += f'<div class="stat-box"><strong>{total_groups}</strong><br>Total Haplotypes</div>\n'
    html += f'<div class="stat-box warning"><strong>{multi_species}</strong><br>Multi-Species Haplotypes</div>\n'
    html += f'<div class="stat-box alert"><strong>{ambiguous}</strong><br>Ambiguous Haplotypes (&lt;70% primary)</div>\n'
    html += '</div>\n'

    # Filter options
    html += '<div class="filter-container">\n'
    html += '<label for="comp-filter">Filter:</label>\n'
    html += '<select id="comp-filter" onchange="filterCompositionTable()">\n'
    html += '<option value="all">All Haplotypes</option>\n'
    html += '<option value="multi">Multi-Species Only</option>\n'
    html += '<option value="ambiguous">Ambiguous Only</option>\n'
    html += '</select>\n'
    html += '</div>\n'

    # Download note
    html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Full table: <code>haplotypes.csv</code> (<code>labels</code>, <code>majority_label</code>, <code>majority_fraction</code>)</p>\n'

    # Interactive table
    html += '<div class="table-container">\n'
    html += '<table id="composition-table" class="data-table">\n'
    html += '<thead><tr>\n'
    html += '<th>Haplotype</th>\n'
    html += '<th>Total Samples</th>\n'
    html += '<th># Species</th>\n'
    html += '<th>Primary Species</th>\n'
    html += '<th>Primary %</th>\n'
    html += '<th>Species Composition</th>\n'
    html += '<th>Flags</th>\n'
    html += '</tr></thead>\n'
    html += '<tbody>\n'

    for _, row in comp_df.iterrows():
        # Apply row highlighting for multi-species/ambiguous
        row_class = ""
        if row['is_ambiguous']:
            row_class = ' class="alert-row"'
        elif row['is_multi_species']:
            row_class = ' class="warning-row"'

        flags = []
        if row['is_multi_species']:
            flags.append('<span class="badge badge-warning">Multi-species</span>')
        if row['is_ambiguous']:
            flags.append('<span class="badge badge-alert">Ambiguous</span>')
        flags_html = " ".join(flags) if flags else "—"

        html += f'<tr{row_class} data-multi="{row["is_multi_species"]}" data-ambiguous="{row["is_ambiguous"]}">\n'
        html += f'<td><strong>{row["haplotype_id"]}</strong></td>\n'
        html += f'<td>{row["n_total_samples"]}</td>\n'
        html += f'<td>{row["n_species"]}</td>\n'
        html += f'<td><em>{row["primary_species"]}</em></td>\n'
        html += f'<td>{row["primary_species_pct"]}</td>\n'
        html += f'<td class="composition-cell">{row["species_composition"]}</td>\n'
        html += f'<td>{flags_html}</td>\n'
        html += '</tr>\n'

    html += '</tbody>\n'
    html += '</table>\n'
    html += '</div>\n'

    # Legend
    html += '<div style="margin-top: 20px; padding: 15px; background: #f8f9fa; border-radius: 4px;">\n'
    html += '<h4 style="margin-top: 0;">Understanding the Flags:</h4>\n'
    html += '<ul style="margin-bottom: 0;">\n'
    html += '<li><strong>Multi-species:</strong> This consensus group contains sequences reported under multiple species names (each >1% abundance)</li>\n'
    html += '<li><strong>Ambiguous:</strong> No single species accounts for >70% of the sequences in this consensus group</li>\n'
    html += '</ul>\n'
    html += '<p style="margin-top: 10px; margin-bottom: 0; font-size: 0.9em; color: #666;">'
    html += '<strong>Note:</strong> These counts describe the records\u2019 BOLD labels. Label conflicts within haplotypes, '
    html += 'with the evidence for each record, are listed in <code>label_conflicts.csv</code>.'
    html += '</p>\n'
    html += '</div>\n'

    html += '</div>\n'  # Close section
    return html


def _build_quality_control_section(quality_dir: Path) -> str:
    """Build quality control summary section."""
    html = '<div class="section" id="section-quality-control">\n'
    html += '<h2>Quality Control</h2>\n'

    if not quality_dir.exists():
        html += '<p class="alert alert-info">Quality control outputs not found.</p>\n'
        html += '</div>\n'
        return html

    summary_df = pd.DataFrame()
    conflicts_df = pd.DataFrame()
    hits_df = pd.DataFrame()

    summary_file = quality_dir / "haplotype_label_summary.csv"
    conflicts_file = quality_dir / "label_conflicts.csv"
    hits_file = quality_dir / "depositor_keyword_hits.csv"

    try:
        if summary_file.exists():
            summary_df = pd.read_csv(summary_file)
        if conflicts_file.exists():
            conflicts_df = pd.read_csv(conflicts_file)
        if hits_file.exists():
            hits_df = pd.read_csv(hits_file)
    except Exception as e:
        html += f'<p class="alert alert-warning">Could not load label outputs: {e}</p>\n'
        html += '</div>\n'
        return html

    total_groups = len(summary_df)
    with_conflict = int((summary_df['n_label_conflicts'] > 0).sum()) if 'n_label_conflicts' in summary_df.columns else 0
    n_missing = int(summary_df['n_missing_label'].sum()) if 'n_missing_label' in summary_df.columns else 0

    html += '<p>Counts compare each haplotype\u2019s member records by their BOLD labels; they do not decide which label is correct. Definitions: <code>README.md</code>.</p>\n'
    html += '<div class="metric-grid">\n'
    html += '  <div class="metric-card"><div class="title">Haplotypes</div>'
    html += f'  <div class="value">{_format_number(total_groups, 0)}</div>'
    html += '  <div class="description">With member labels summarized</div></div>\n'

    html += '  <div class="metric-card"><div class="title">With a Label Conflict</div>'
    html += f'  <div class="value">{_format_number(with_conflict, 0)}</div>'
    pct = (with_conflict / total_groups * 100) if total_groups else 0
    html += f'  <div class="description">{_format_percentage(pct)} of haplotypes</div></div>\n'

    html += '  <div class="metric-card"><div class="title">Label-Conflict Records</div>'
    html += f'  <div class="value">{_format_number(len(conflicts_df), 0)}</div>'
    html += '  <div class="description">See label_conflicts.csv</div></div>\n'

    html += '  <div class="metric-card"><div class="title">Members Without a Label</div>'
    html += f'  <div class="value">{_format_number(n_missing, 0)}</div>'
    html += '  <div class="description">Not counted in majority fractions</div></div>\n'

    html += '  <div class="metric-card"><div class="title">Depositor-Note Keyword Hits</div>'
    html += f'  <div class="value">{_format_number(len(hits_df), 0)}</div>'
    html += '  <div class="description">Notes to read; not interpreted</div></div>\n'
    html += '</div>\n'

    if not summary_df.empty:
        html += '<h3>Haplotypes with Label Conflicts</h3>\n'
        html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Full table: <code>haplotypes.csv</code></p>\n'
        flagged = summary_df[summary_df['n_label_conflicts'] > 0] if 'n_label_conflicts' in summary_df.columns else summary_df
        keep_cols = [c for c in ['haplotype_id', 'n_members', 'n_labeled', 'n_missing_label', 'majority_label',
                                 'majority_fraction', 'labels', 'n_label_conflicts'] if c in flagged.columns]
        if not flagged.empty and keep_cols:
            html += _dataframe_to_html(flagged[keep_cols], max_rows=20)
        elif flagged.empty:
            html += '<p class="alert alert-success">No haplotype has members with conflicting labels.</p>\n'

    if not conflicts_df.empty:
        html += '<h3>Label-Conflict Records</h3>\n'
        html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Details: <code>label_conflicts.csv</code></p>\n'
        keep_cols = [c for c in ['processid', 'label', 'haplotype_id', 'majority_label', 'n_majority', 'n_labeled',
                                 'majority_fraction', 'nearest_haplotype_with_label',
                                 'distance_to_nearest_with_label'] if c in conflicts_df.columns]
        if keep_cols:
            html += _dataframe_to_html(conflicts_df[keep_cols], max_rows=15)

    if not hits_df.empty:
        html += '<h3>Depositor-Note Keyword Hits</h3>\n'
        html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Details: <code>records.csv</code> (<code>depositor_keyword_hits</code>); with --extended, <code>extended/labels/depositor_keyword_hits.csv</code></p>\n'
        keep_cols = [c for c in ['processid', 'haplotype_id', 'label', 'keyword_hits', 'depositor_note'] if c in hits_df.columns]
        if keep_cols:
            html += _dataframe_to_html(hits_df[keep_cols], max_rows=15)

    html += '</div>\n'
    return html


def _load_provenance(run_dir: Path, provenance: Optional[dict]) -> dict:
    """The run's provenance: the dict passed by the pipeline, else run_provenance.json."""
    if provenance is not None:
        return provenance
    try:
        return json.loads((Path(run_dir) / "run_provenance.json").read_text())
    except (OSError, ValueError):
        return {}


def _build_parameters_section(output_dir: Path, organism: str,
                              provenance: Optional[dict] = None) -> str:
    """Build pipeline parameters section."""
    html = '<div class="section" id="section-parameters">\n'
    html += '<h2>Pipeline Parameters</h2>\n'

    # Parameters come from the run's provenance (the fully resolved configuration)
    prov = _load_provenance(output_dir, provenance)
    if not prov:
        html += '<p class="alert alert-info">Parameter information not available</p>\n'
        html += '</div>\n'
        return html

    try:
        cfg = prov.get('config', {})
        applied = prov.get('applied', {})

        html += '<p>The following parameters were used for this analysis:</p>\n'
        html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 10px;">Full configuration, tool versions and input checksum: <code>run_provenance.json</code></p>\n'
        html += '<table class="data-table">\n'
        html += '<thead><tr><th>Parameter</th><th>Value</th><th>Description</th></tr></thead>\n'
        html += '<tbody>\n'

        def _get_nested(d, *keys, default='N/A'):
            for k in keys:
                if isinstance(d, dict):
                    d = d.get(k, default)
                else:
                    return default
            return d

        def _row(label, value, description, pct=None, unit=''):
            html_row = f'<tr><td><strong>{label}</strong></td><td><code>{value}</code>{unit}'
            if pct and isinstance(value, (int, float)) and not isinstance(value, bool):
                html_row += f' ({value*100:{pct}}%)'
            return html_row + f'</td><td>{description}</td></tr>\n'

        html += _row('Singleton error filter', _get_nested(cfg, 'haplotype', 'min_singleton_distance'),
                     'Singleton haplotypes within this p-distance of another haplotype are treated as '
                     'sequencing or PCR errors (--singleton-distance).', pct='.1f')
        html += _row('Error singletons assigned to nearest haplotype',
                     _get_nested(cfg, 'haplotype', 'assign_corrected'),
                     'Records of error singletons get status assigned_corrected (unless equidistant).')
        html += _row('Suspect singleton threshold', _get_nested(cfg, 'haplotype', 'max_singleton_distance'),
                     'Retained singletons farther than this from their nearest haplotype are flagged suspect.',
                     pct='.1f')
        html += _row('Core region column coverage', _get_nested(cfg, 'core_region', 'core_min_coverage'),
                     'A position is part of the core region if at least this fraction of sequences covers it.',
                     pct='.0f')
        html += _row('Per-sequence core cutoff', _get_nested(cfg, 'core_region', 'min_sequence_coverage'),
                     'A sequence needs at least this fraction of the reference core length '
                     f'({applied.get("min_sequence_length_bp", "N/A"):.0f} bp applied).'
                     if isinstance(applied.get('min_sequence_length_bp'), (int, float))
                     else 'A sequence needs at least this fraction of the reference core length.',
                     pct='.0f')
        html += _row('Minimum sequence length', _get_nested(cfg, 'qc', 'min_raw_length_abs'),
                     'Absolute minimum raw sequence length for QC filtering.', unit=' bp')
        html += _row('Maximum N fraction', _get_nested(cfg, 'qc', 'max_raw_N_fraction'),
                     'Maximum fraction of ambiguous bases (N) allowed.', pct='.1f')
        html += _row('Threads', cfg.get('n_threads', 'N/A'), 'Number of parallel processing threads used.')
        html += _row('Phylogenetic tree', _get_nested(cfg, 'phylogenetic', 'build_tree'),
                     'Whether a phylogenetic tree was built.')

        html += '</tbody>\n'
        html += '</table>\n'

    except Exception as e:
        html += f'<p class="alert alert-warning">Could not load parameters: {e}</p>\n'

    html += '</div>\n'
    return html


def _build_geographic_section(annotated_df: pd.DataFrame) -> str:
    """Build geographic distribution section."""
    html = '<div class="section" id="section-geographic">\n'
    html += '<h2>Geographic Distribution</h2>\n'
    html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 15px;">Coordinates, coordinate status and regions per record: <code>records.csv</code></p>\n'

    if 'ocean_basin' not in annotated_df.columns:
        html += '<p class="alert alert-info">No geographic data available</p>\n'
        html += '</div>\n'
        return html

    # Count total samples and samples with Unknown geography
    total_samples = len(annotated_df)
    unknown_geography = (annotated_df['ocean_basin'].fillna('Unknown') == 'Unknown').sum()
    samples_with_geography = total_samples - unknown_geography

    # Add summary information
    html += '<div class="alert alert-info">\n'
    html += f'<strong>Geographic Analysis Summary:</strong> {_format_number(samples_with_geography, 0)} of '
    html += f'{_format_number(total_samples, 0)} samples '
    html += f'({_format_percentage(samples_with_geography / total_samples * 100)}) have defined geographic locations. '
    if unknown_geography > 0:
        html += f'{_format_number(unknown_geography, 0)} samples '
        html += f'({_format_percentage(unknown_geography / total_samples * 100)}) '
        html += 'with unknown or missing geography were excluded from geographic analyses.\n'
    html += '</div>\n'

    # Region distribution (excluding Unknown)
    html += '<h3>Sample Distribution by Geographic Region</h3>\n'

    # Filter out Unknown geography for the table
    basin_series = annotated_df['ocean_basin'].fillna('Unknown')
    basin_counts = basin_series[basin_series != 'Unknown'].value_counts().reset_index()

    if len(basin_counts) > 0:
        basin_counts.columns = ['Geographic Region', 'Sample Count']
        # Calculate percentages based on samples with defined geography
        basin_counts['Percentage'] = (basin_counts['Sample Count'] / samples_with_geography * 100).round(1)
        html += _dataframe_to_html(basin_counts, max_rows=20)
    else:
        html += '<p class="alert alert-warning">No samples with defined geographic region assignments</p>\n'

    # Haplotype × Region distribution (excluding Unknown)
    if 'haplotype_id' in annotated_df.columns:
        html += '<h3>Haplotypes per Geographic Region</h3>\n'

        # Filter out Unknown geography before creating crosstab
        df_with_geography = annotated_df[annotated_df['ocean_basin'].fillna('Unknown') != 'Unknown'].copy()

        if len(df_with_geography) > 0:
            # Create crosstab
            ct = pd.crosstab(
                df_with_geography['haplotype_id'].fillna('Unassigned'),
                df_with_geography['ocean_basin']
            )

            # Convert to long format for better display
            ct_reset = ct.reset_index()
            html += _dataframe_to_html(ct_reset, max_rows=30)
        else:
            html += '<p class="alert alert-warning">No genotyped samples with defined geographic region assignments</p>\n'

    # Missing geography
    if 'lat' in annotated_df.columns and 'lon' in annotated_df.columns:
        missing_coords = annotated_df[annotated_df[['lat', 'lon']].isna().any(axis=1)]
        n_missing = len(missing_coords)
        pct_missing = (n_missing / len(annotated_df) * 100)

        if n_missing > 0:
            html += '<div class="alert alert-warning">\n'
            html += f'<strong>Missing Geography:</strong> {_format_number(n_missing, 0)} samples '
            html += f'({_format_percentage(pct_missing)}) do not have valid coordinate data.\n'
            html += '</div>\n'

    html += '</div>\n'  # Close section
    return html


def _build_metadata_summary_section(
    output_dir: Path,
    organism: str,
    resolve=None,
) -> str:
    """
    Build metadata summarization section for HTML report.

    Includes coverage statistics, per-field categorical tabulations, and
    temporal trends. Hypothesis testing (chi-square) was removed; see
    metadata_summary.py module docstring for rationale.

    Parameters
    ----------
    output_dir : Path
        Base output directory
    organism : str
        Organism name

    Returns
    -------
    str
        HTML content for metadata summary section
    """
    html = '<div class="section" id="section-metadata-summary">\n'
    html += '<h2>Metadata Summary</h2>\n'
    resolve = resolve or _relative_link(output_dir)

    metadata_dir = output_dir / 'metadata_summary'
    viz_dir = output_dir / 'visualization' / 'metadata'
    viz_png_dir = viz_dir / 'png'

    def _find_viz_image(basename):
        """Find a visualization image, checking png/ subdir first, then base dir."""
        for d in [viz_png_dir, viz_dir / 'pdf', viz_dir]:
            p = d / basename
            if p.exists():
                return p
        return viz_dir / basename  # Default path even if not found

    # Check if metadata summarization was performed
    if not metadata_dir.exists():
        html += '<p class="alert alert-info">Metadata summarization not performed. '
        html += 'It runs by default; use <code>--no-metadata-summary</code> to disable.</p>\n'
        html += '</div>\n'
        return html

    html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 15px;">'
    html += 'Field coverage: <code>metadata_coverage.csv</code>; with --extended, cross-tabulations and collection-date summaries in <code>extended/metadata/</code></p>\n'

    # 1. Coverage Overview
    coverage_json = metadata_dir / f'{organism}_metadata_coverage.json'
    if coverage_json.exists():
        html += '<h3>Metadata Field Coverage</h3>\n'
        html += '<p>Coverage indicates the percentage of samples with non-null values for each metadata field.</p>\n'

        try:
            with open(coverage_json, 'r') as f:
                coverage_data = json.load(f)

            total_samples = coverage_data.get('total_samples', 0)
            html += f'<p><strong>Total Samples:</strong> {_format_number(total_samples, 0)}</p>\n'

            # Create coverage table
            fields_data = coverage_data.get('fields', {})
            if fields_data:
                rows = []
                for field, stats in fields_data.items():
                    rows.append({
                        'Field': field,
                        'Coverage': f"{stats.get('pct_coverage', 0):.1f}%",
                        'Samples': _format_number(stats.get('n_with_value', 0), 0),
                        'Unique Values': stats.get('n_unique', 0),
                        'Top Values': ', '.join(stats.get('top_values', [])[:3])
                    })
                    # Add warning if present
                    if stats.get('warning'):
                        rows[-1]['Notes'] = stats['warning']

                coverage_df = pd.DataFrame(rows)
                html += _dataframe_to_html(coverage_df, max_rows=20)

                # Check for data quality warnings
                warnings = [f"{f}: {s['warning']}" for f, s in fields_data.items() if s.get('warning')]
                if warnings:
                    html += '<div class="alert alert-warning">\n'
                    html += '<strong>Data Quality Notes:</strong>\n<ul>\n'
                    for w in warnings:
                        html += f'<li>{w}</li>\n'
                    html += '</ul>\n</div>\n'

        except Exception as e:
            html += f'<p class="alert alert-warning">Could not load coverage data: {e}</p>\n'

        # Embed coverage visualization if available
        coverage_img = _find_viz_image(f'{organism}_metadata_coverage.png')
        if coverage_img.exists():
            src = resolve(coverage_img)
            if src:
                html += '<div class="viz-container">\n'
                html += f'<img src="{src}" alt="Metadata Coverage" class="viz-image" style="max-width: 800px;">\n'
                html += '</div>\n'

    # 2. Collection dates (--extended)
    temporal_summary = metadata_dir / f'{organism}_temporal_summary.csv'
    first_sampled_csv = metadata_dir / f'{organism}_haplotype_first_sampled.csv'

    if temporal_summary.exists() or first_sampled_csv.exists():
        html += '<h3>Collection Dates</h3>\n'
        html += ('<p>Collection dates show when and where specimens were sampled and deposited '
                 'in BOLD, not when haplotypes arose or changed in frequency.</p>\n')

        if first_sampled_csv.exists():
            html += '<h4>First Sampled</h4>\n'
            try:
                first_df = pd.read_csv(first_sampled_csv)
                if not first_df.empty:
                    display_df = first_df[['haplotype_sp', 'first_sampled', 'first_location', 'n_total_samples']].copy()
                    display_df.columns = ['Haplotype', 'First Sampled', 'Location', 'Total Samples']
                    html += _dataframe_to_html(display_df, max_rows=15)
            except Exception as e:
                html += f'<p class="alert alert-warning">Could not load first-sampled data: {e}</p>\n'

        temporal_img = _find_viz_image(f'{organism}_temporal_composition.png')
        if temporal_img.exists():
            src = resolve(temporal_img)
            if src:
                html += '<h4>Records by Collection Year</h4>\n'
                html += '<div class="viz-container">\n'
                html += f'<img src="{src}" alt="Records by collection year" class="viz-image" style="max-width: 800px;">\n'
                html += '</div>\n'

        timeline_plot = _find_viz_image(f'{organism}_haplotype_sampling_timeline.png')
        if timeline_plot.exists():
            src = resolve(timeline_plot)
            if src:
                html += '<h4>Sampling Timeline</h4>\n'
                html += '<div class="viz-container">\n'
                html += f'<img src="{src}" alt="Sampling timeline" class="viz-image" style="max-width: 800px;">\n'
                html += '</div>\n'

    # 4. Metadata Heatmap
    heatmap_img = _find_viz_image(f'{organism}_metadata_heatmap.png')
    if heatmap_img.exists():
        html += '<h3>Metadata Coverage Heatmap</h3>\n'
        html += '<p>Coverage of metadata fields across haplotypes.</p>\n'
        src = resolve(heatmap_img)
        if src:
            html += '<div class="viz-container">\n'
            html += f'<img src="{src}" alt="Metadata Heatmap" class="viz-image" style="max-width: 100%;">\n'
            html += '</div>\n'

    # 5. Selected Categorical Analyses
    categorical_files = list(metadata_dir.glob(f'{organism}_*_tabulation.csv'))
    categorical_files = [f for f in categorical_files if 'association' not in f.name and 'temporal' not in f.name]

    if categorical_files:
        html += '<h3>Categorical Field Analyses</h3>\n'
        html += '<p>Haplotype distribution across selected metadata categories. '
        html += 'Full tabulations in <code>extended/metadata/</code>.</p>\n'

        # Create tabs for each categorical analysis
        html += '<div class="subtabs">\n'
        for idx, cat_file in enumerate(categorical_files[:6]):  # Limit to 6 tabs
            field_name = cat_file.stem.replace(f'{organism}_', '').replace('_tabulation', '')
            active = ' active' if idx == 0 else ''
            tab_id = f'meta-cat-{idx}'
            html += f'<button class="subtab-button{active}" onclick="switchSubtab(\'{tab_id}\')">{field_name}</button>\n'
        html += '</div>\n'

        for idx, cat_file in enumerate(categorical_files[:6]):
            field_name = cat_file.stem.replace(f'{organism}_', '').replace('_tabulation', '')
            active = ' active' if idx == 0 else ''
            tab_id = f'meta-cat-{idx}'

            html += f'<div id="{tab_id}" class="subtab-content{active}">\n'

            try:
                cat_df = pd.read_csv(cat_file)
                if not cat_df.empty:
                    # Show top entries
                    display_df = cat_df.head(20).copy()
                    if 'pct_of_haplotype' in display_df.columns:
                        display_df['pct_of_haplotype'] = display_df['pct_of_haplotype'].apply(lambda x: f'{x:.1f}%')
                    if 'pct_of_value' in display_df.columns:
                        display_df['pct_of_value'] = display_df['pct_of_value'].apply(lambda x: f'{x:.1f}%')
                    html += _dataframe_to_html(display_df, max_rows=20)

                    # Show corresponding visualization if available
                    field_clean = field_name.replace('/', '_')
                    viz_file = _find_viz_image(f'{organism}_{field_clean}_by_haplotype.png')
                    if viz_file.exists():
                        src = resolve(viz_file)
                        if src:
                            html += '<div class="viz-container" style="margin-top: 15px;">\n'
                            html += f'<img src="{src}" alt="{field_name} Distribution" class="viz-image" style="max-width: 100%;">\n'
                            html += '</div>\n'
            except Exception as e:
                html += f'<p class="alert alert-warning">Could not load {field_name} analysis: {e}</p>\n'

            html += '</div>\n'

    html += '</div>\n'  # Close section
    return html


def _relative_link(base_dir: Path):
    """Link resolver for a report written in base_dir: path relative to it."""
    import os
    base = Path(base_dir).resolve()
    return lambda p: os.path.relpath(Path(p).resolve(), base).replace(os.sep, "/")


def _link_resolver(links: dict):
    """Link resolver for staged files: their path in the run folder, or None
    when the run does not keep the file (P2: the report links, it does not embed)."""
    table = {Path(src).resolve(): Path(rel).as_posix() for src, rel in links.items()}
    return lambda p: table.get(Path(p).resolve())


def _encode_image_to_base64(image_path: Path) -> Optional[str]:
    """
    Encode an image file to base64 for embedding in HTML.

    Parameters
    ----------
    image_path : Path
        Path to image file

    Returns
    -------
    Optional[str]
        Base64-encoded image data URI, or None if encoding failed
    """
    import base64

    if not image_path.exists():
        return None

    try:
        with open(image_path, 'rb') as f:
            image_data = f.read()

        # Determine MIME type from extension
        ext = image_path.suffix.lower()
        mime_types = {
            '.png': 'image/png',
            '.jpg': 'image/jpeg',
            '.jpeg': 'image/jpeg',
            '.gif': 'image/gif',
            '.svg': 'image/svg+xml'
        }
        mime_type = mime_types.get(ext, 'image/png')

        # Encode to base64
        encoded = base64.b64encode(image_data).decode('utf-8')
        return f"data:{mime_type};base64,{encoded}"

    except Exception as e:
        logger.warning(f"Failed to encode image {image_path}: {e}", extra={"step_failure": "12 HTML report"})
        return None


def _build_visualizations_section(
    output_dir: Path,
    organism: str,
    resolve=None,
) -> str:
    """
    Build visualizations section with embedded images.

    Parameters
    ----------
    output_dir : Path
        Base output directory
    organism : str
        Organism name

    Returns
    -------
    str
        HTML content for visualizations section
    """
    html = '<div class="section" id="section-visualizations">\n'
    html += '<h2>Visualizations</h2>\n'
    resolve = resolve or _relative_link(output_dir)
    html += '<p style="color: #666; font-size: 0.9em; margin-bottom: 15px;">Figures are linked from the run folder (<code>geography/</code>, <code>phylogeny/</code>, <code>extended/</code>); open this report from that folder.</p>\n'

    # Add global color picker container
    html += '<div id="global-color-picker" class="plot-controls" style="margin-bottom: 20px; max-width: 600px;">\n'
    html += '<!-- Global color picker will be populated by JavaScript -->\n'
    html += '</div>\n'

    viz_dir = output_dir / 'visualization'
    viz_pdf_dir = viz_dir / 'pdf'
    viz_svg_dir = viz_dir / 'svg'
    viz_json_dir = viz_dir / 'json'

    if not viz_dir.exists():
        html += '<p class="alert alert-info">No visualizations available</p>\n'
        html += '</div>\n'
        return html

    # Define visualization categories and their display names (in display order)
    viz_categories = [
        {
            'title': 'Phylogenetic Tree',
            'pattern': f'{organism}_tree.png',
            'directory': 'phylogenetic',
            'description': 'Phylogenetic tree showing relationships between haplotypes',
            'note': f'Newick tree files: <code>phylogeny/tree.nwk</code> (haplotype IDs) and <code>phylogeny/tree_labeled.nwk</code> (haplotype labels), for opening in tree editors such as <a href="https://treeviewer.org/" target="_blank">TreeViewer</a> (Bianchini &amp; Sánchez-Baracaldo, 2024) for re-rooting and customization'
        },
        {
            'title': 'Relative Abundance by Geographic Region',
            'pattern': f'{organism}_distribution_bar.png',
            'json_pattern': f'{organism}_distribution_bar_data.json',
            'description': 'Relative abundance of haplotypes across geographic regions'
        },
        {
            'title': 'Total Abundance by Geographic Region',
            'pattern': f'{organism}_totaldistribution_bar.png',
            'json_pattern': f'{organism}_totaldistribution_bar_data.json',
            'description': 'Total sample counts of haplotypes across geographic regions'
        },
        {
            'title': 'Haplotype Abundance by Region',
            'pattern': f'{organism}_distribution_bar_faceted.png',
            'json_pattern': f'{organism}_distribution_bar_faceted_data.json',
            'description': 'Total sample counts faceted by haplotype across geographic regions'
        },
        {
            'title': 'Species Abundance by Region',
            'pattern': f'{organism}_distribution_bar_species_faceted.png',
            'json_pattern': f'{organism}_distribution_bar_species_faceted_data.json',
            'description': 'Total sample counts faceted by assigned species with haplotypes stacked across geographic regions'
        },
        {
            'title': 'Distribution Map',
            'pattern': f'{organism}_distribution_map.png',
            'json_pattern': f'{organism}_distribution_map_data.json',
            'description': 'Geographic distribution of samples'
        },
        {
            'title': 'Distribution Map (Species-Faceted)',
            'pattern': f'{organism}_distribution_map_species_faceted.png',
            'description': 'Geographic distribution faceted by assigned species'
        },
        {
            'title': 'Distribution Map (Haplotype-Faceted)',
            'pattern': f'{organism}_distribution_map_faceted.png',
            'description': 'Geographic distribution faceted by haplotype'
        }
    ]

    # Filter to available visualizations
    available_viz = []
    for viz in viz_categories:
        # Determine which directory to look in
        base_search_dir = output_dir / viz.get('directory', 'visualization')

        # Allow SVG/PNG fallbacks (ESV workflow exports SVG)
        # Check both new structure (pdf/svg subdirs) and legacy (flat structure)
        # Note: PDFs are excluded because they cannot be embedded in HTML <img> tags
        base_pattern = viz['pattern']
        base_name = Path(base_pattern).stem
        candidate_paths = [
            # New structure (PNG and SVG only — embeddable in HTML)
            viz_pdf_dir / Path(base_pattern).with_suffix('.png').name,
            viz_svg_dir / Path(base_pattern).with_suffix('.svg').name,
            # Legacy flat structure
            viz_dir / base_pattern,
            viz_dir / Path(base_pattern).with_suffix('.svg'),
            viz_dir / Path(base_pattern).with_suffix('.png'),
            # Non-visualization directories (e.g., phylogenetic)
            base_search_dir / base_pattern,
            base_search_dir / Path(base_pattern).with_suffix('.svg'),
            base_search_dir / Path(base_pattern).with_suffix('.png'),
        ]
        # Link to the first candidate that the run folder keeps (not embedded)
        image_path, image_src = next(((p, resolve(p)) for p in candidate_paths
                                      if p.exists() and resolve(p)), (None, None))

        if image_path:
            files = []
            for ext in ('.pdf', '.svg'):
                name = Path(base_pattern).with_suffix(ext).name
                staged = next((d / name for d in (viz_pdf_dir, viz_svg_dir, viz_dir, base_search_dir)
                               if (d / name).exists() and resolve(d / name)), None)
                if staged is not None:
                    files.append(resolve(staged))
            if image_src:
                viz_data = {**viz, 'image_src': image_src, 'files': files}

                # Load JSON data if available (for interactive plots)
                if 'json_pattern' in viz:
                    # Check both new and legacy locations
                    json_candidates = [
                        viz_json_dir / viz['json_pattern'],
                        viz_dir / viz['json_pattern'],
                    ]
                    json_path = next((p for p in json_candidates if p.exists()), None)
                    if json_path:
                        try:
                            with open(json_path, 'r') as f:
                                plot_data = json.load(f)
                            viz_data['plot_data'] = plot_data
                            logger.debug(f"Loaded plot data for {viz['title']}: {json_path}")
                        except Exception as e:
                            logger.warning(f"Failed to load plot data from {json_path}: {e}", extra={"step_failure": "12 HTML report"})

                available_viz.append(viz_data)

    if not available_viz:
        html += '<p class="alert alert-info">No visualization images found</p>\n'
        html += '</div>\n'
        return html

    # Create subtabs for visualizations
    html += '<div class="subtabs">\n'
    for idx, viz in enumerate(available_viz):
        active_class = ' active' if idx == 0 else ''
        tab_id = f"viz-tab-{idx}"
        html += f'<button class="subtab-button{active_class}" onclick="switchSubtab(\'{tab_id}\')">{viz["title"]}</button>\n'
    html += '</div>\n'

    # Create subtab content
    for idx, viz in enumerate(available_viz):
        active_class = ' active' if idx == 0 else ''
        tab_id = f"viz-tab-{idx}"
        html += f'<div id="{tab_id}" class="subtab-content{active_class}">\n'
        html += '<div class="viz-container">\n'
        html += f'<h3>{viz["title"]}</h3>\n'
        html += f'<p class="viz-description">{viz["description"]}</p>\n'
        links = ', '.join(f'<a href="{f}"><code>{f}</code></a>' for f in viz['files'])
        html += f'<p style="color: #666; font-size: 0.85em; margin-bottom: 10px;">Files: {links}</p>\n'
        # Add optional note if present
        if 'note' in viz:
            html += f'<p style="color: #2563eb; font-size: 0.85em; margin-bottom: 10px; padding: 8px; background-color: #eff6ff; border-left: 3px solid #2563eb; border-radius: 3px;">{viz["note"]}</p>\n'

        # Add interactive controls for plots with JSON data (except faceted plots)
        if 'plot_data' in viz and viz['plot_data'].get('plot_type') != 'faceted_abundance':
            data_id = f"plot-data-{idx}"
            plot_id = f"plot-{idx}"
            controls_id = f"controls-{idx}"

            # Toggle button between static and interactive
            html += '<div class="plot-mode-toggle">\n'
            html += f'<button class="toggle-btn active" onclick="showStaticPlot(\'{plot_id}\', {idx})">Static Image</button>\n'
            html += f'<button class="toggle-btn" onclick="showInteractivePlot(\'{plot_id}\', {idx})">Interactive Plot</button>\n'
            html += '</div>\n'

            # Static plot (shown by default)
            html += f'<div id="{plot_id}-static" class="plot-view active">\n'
            html += f'<img src="{viz["image_src"]}" alt="{viz["title"]}" class="viz-image">\n'
            html += '</div>\n'

            # Interactive plot container (hidden by default)
            html += f'<div id="{plot_id}-interactive" class="plot-view" style="display: none;">\n'

            # Filter controls
            html += f'<div id="{controls_id}" class="plot-controls">\n'

            # Summary Statistics Panel
            html += f'<div id="stats-panel-{idx}" class="stats-panel">\n'
            html += '<div class="stats-grid">\n'
            html += '<div class="stat-item">\n'
            html += '<div class="stat-label">Genotypes Shown</div>\n'
            html += f'<div class="stat-value" id="stat-genotypes-{idx}">--</div>\n'
            html += '</div>\n'
            html += '<div class="stat-item">\n'
            html += '<div class="stat-label">Total Samples</div>\n'
            html += f'<div class="stat-value" id="stat-samples-{idx}">--</div>\n'
            html += '</div>\n'
            html += '<div class="stat-item">\n'
            html += '<div class="stat-label">% of Total Data</div>\n'
            html += f'<div class="stat-value" id="stat-percentage-{idx}">--</div>\n'
            html += '</div>\n'
            html += '</div>\n'
            html += '</div>\n'

            html += '<div class="controls-section">\n'
            html += '<h4>Filter Options</h4>\n'

            # Threshold slider
            html += '<div class="control-group">\n'
            html += '<label for="threshold-{0}">Minimum Sample Count: <span id="threshold-value-{0}">0</span></label>\n'.format(idx)
            html += '<input type="range" id="threshold-{0}" class="threshold-slider" min="0" max="50" value="0" \n'.format(idx)
            html += f'oninput="updateThreshold({idx}); updatePlot({idx})">\n'
            html += '</div>\n'

            # Wrap basin filter and genotype checkboxes in flex container
            html += '<div class="controls-flex-container">\n'

            # Ocean basin filter (only for distribution maps with ocean basin data)
            if viz['plot_data'].get('plot_type') == 'distribution_map' and viz['plot_data'].get('ocean_basins'):
                html += '<div class="control-group">\n'
                html += '<label>Filter by Ocean Basin:</label>\n'
                html += f'<div class="checkbox-controls">\n'
                html += f'<button class="btn-small" onclick="selectAllBasins({idx})">Select All</button>\n'
                html += f'<button class="btn-small" onclick="deselectAllBasins({idx})">Deselect All</button>\n'
                html += '</div>\n'
                html += f'<div id="basin-checkboxes-{idx}" class="genotype-checkboxes">\n'
                html += '<!-- Basin checkboxes will be populated by JavaScript -->\n'
                html += '</div>\n'
                html += '</div>\n'

            # Genotype checkboxes
            html += '<div class="control-group">\n'
            html += '<label>Select Genotypes:</label>\n'
            html += f'<div class="checkbox-controls">\n'
            html += f'<button class="btn-small" onclick="selectAllGenotypes({idx})">Select All</button>\n'
            html += f'<button class="btn-small" onclick="deselectAllGenotypes({idx})">Deselect All</button>\n'
            html += '</div>\n'
            html += f'<div id="genotype-checkboxes-{idx}" class="genotype-checkboxes">\n'
            html += '<!-- Checkboxes will be populated by JavaScript -->\n'
            html += '</div>\n'
            html += '</div>\n'

            html += '</div>\n'  # Close controls-flex-container

            # Download options
            html += '<div class="control-group">\n'
            html += '<label>Export Options:</label>\n'
            html += '<div class="download-options">\n'
            html += f'<button class="btn-download" onclick="downloadPlot({idx}, \'png\')">📥 PNG</button>\n'
            html += f'<button class="btn-download" onclick="downloadPlot({idx}, \'svg\')">📥 SVG</button>\n'
            html += f'<button class="btn-download" onclick="downloadCSV({idx})">📥 CSV Data</button>\n'
            html += '</div>\n'
            html += '<div class="resolution-control">\n'
            html += '<label for="resolution-{0}" class="small-label">Image Resolution:</label>\n'.format(idx)
            html += '<select id="resolution-{0}" class="resolution-select">\n'.format(idx)
            html += '<option value="800x600">Standard (800×600)</option>\n'
            html += '<option value="1200x900" selected>High (1200×900)</option>\n'
            html += '<option value="1600x1200">Print (1600×1200)</option>\n'
            html += '<option value="2400x1800">Publication (2400×1800)</option>\n'
            html += '</select>\n'
            html += '</div>\n'
            html += '</div>\n'

            html += '</div>\n'  # controls-section
            html += '</div>\n'  # plot-controls

            # Plotly chart container
            html += f'<div id="{plot_id}" class="plotly-chart"></div>\n'

            html += '</div>\n'  # plot-view interactive

            # Embed plot data as JSON
            json_str = json.dumps(viz['plot_data'])
            html += f'<script type="application/json" id="{data_id}">\n'
            html += json_str
            html += '\n</script>\n'
        else:
            # No interactive data or faceted plot - just show static image
            html += f'<img src="{viz["image_src"]}" alt="{viz["title"]}" class="viz-image">\n'

        html += '</div>\n'  # viz-container
        html += '</div>\n'  # subtab-content

    html += '</div>\n'
    return html


def _build_methods_section(output_dir: Path, organism: str, version: str,
                           provenance: Optional[dict] = None) -> str:
    """
    Methods for the haplotype (exact sequence variant) workflow, with the
    values this run used: the pipeline's provenance, else run_provenance.json.
    """
    prov = _load_provenance(output_dir, provenance)
    cfg = prov.get('config', {})
    applied = prov.get('applied', {})
    tools = prov.get('tools', {})
    summary = applied.get('assignment_summary', {})
    qc = cfg.get('qc', {})
    coi = cfg.get('coi_validation', {})
    core = cfg.get('core_region', {})
    hap = cfg.get('haplotype', {})
    phy = cfg.get('phylogenetic', {})
    tax = cfg.get('taxonomy', {})
    notes = prov.get('data_notes', {})

    def v(value, fmt="{}", missing="not recorded"):
        return fmt.format(value) if value is not None else missing

    def row(name, value):
        return f'<tr><td><strong>{name}</strong></td><td>{value}</td></tr>\n'

    def table(rows):
        return '<table class="data-table">\n<tbody>\n' + ''.join(rows) + '</tbody>\n</table>\n'

    code = coi.get('mitochondrial_code')
    code_name = {2: "vertebrate mitochondrial", 5: "invertebrate mitochondrial"}.get(code, "NCBI genetic code")
    tree_built = applied.get('n_tree_tips') is not None
    rooting = applied.get('tree_rooting') or {}
    outgroup = html_escape(str(rooting.get('outgroup', '')))
    kept = '; support values kept on their branches'
    rooting_text = {
        'midpoint': 'Midpoint: the middle of the longest path between two tips' + kept,
        'outgroup_label': f'On the outgroup tip {outgroup}' + kept,
        'outgroup_taxon': f'On the clade of haplotypes assigned to {outgroup}' + kept,
        'outgroup_fasta': f'On the outgroup sequences in {outgroup}' + kept,
        'none': 'Unrooted, as written by FastTree',
    }.get(rooting.get('method'), 'not recorded')
    curation = applied.get('taxon_curation') or {}
    geo = applied.get('geographic_analysis', False)
    geo_cat = applied.get('geo_category') or 'ocean_basin'
    custom_shp = applied.get('custom_shapefile')
    d = hap.get('min_singleton_distance')

    tabs = [('Overview', 'overview'), ('Quality Control', 'qc'), ('Haplotypes', 'haplotypes'),
            ('Record Statuses', 'statuses'), ('Taxonomy & Labels', 'taxonomy'),
            ('Phylogeny', 'phylogenetic'), ('Geography', 'geographic'), ('Software', 'software'),
            ('Methods Statement', 'statement'), ('References', 'references')]

    html = '<div class="section" id="section-methods">\n<h2>Methods</h2>\n'
    html += ('<p class="section-description">What this run did, with the values it used '
             '(from <code>run_provenance.json</code>).</p>\n')
    html += '<div class="subtabs">\n'
    for idx, (title, tab_id) in enumerate(tabs):
        active = ' active' if idx == 0 else ''
        html += f'<button class="subtab-button{active}" onclick="switchSubtab(\'methods-{tab_id}\')">{title}</button>\n'
    html += '</div>\n'

    def tab(tab_id, title, body, active=False):
        return (f'<div id="methods-{tab_id}" class="subtab-content{" active" if active else ""}">\n'
                f'<h3>{title}</h3>\n{body}</div>\n')

    inp = prov.get('input', {})
    html += tab('overview', 'Overview', table([
        row('BOLDGenotyper', f'v{prov.get("boldgenotyper_version", version)}'),
        row('Run started', v(prov.get('started_at'))),
        row('Input', f'<code>{v(inp.get("path"))}</code> ({v(inp.get("n_rows"), "{:,}")} rows; md5 {v(inp.get("md5"))})'),
        row('BOLD download date', v(inp.get('download_date'))),
        row('Command', f'<code>{v(prov.get("command_line"))}</code>'),
    ]), active=True)

    html += tab('qc', 'Quality Control', table([
        row('Duplicate processids', f'{v(notes.get("n_duplicate_processids"), "{:,}")} removed (first record kept)'),
        row('Orientation and reading frame',
            f'Sequences oriented and checked for an open reading frame with genetic code {v(code)} '
            f'({code_name}): ORF coverage at least {v(coi.get("orf_min_coverage"))}, at most '
            f'{v(coi.get("orf_max_internal_stops"))} internal stop codons'),
        row('Length', f'At least {v(qc.get("min_raw_length_abs"))} bp and at least '
                      f'{v(qc.get("min_raw_length_frac_of_median"))} of the median length'),
        row('Ambiguous bases', f'N fraction at most {v(qc.get("max_raw_N_fraction"))}'),
        row('Valid reading frame required', v(qc.get('require_valid_orf'))),
        row('Records passing QC', f'{v(summary.get("n_qc_pass"), "{:,}")} of {v(summary.get("n_unique"), "{:,}")} unique records'),
    ]))

    html += tab('haplotypes', 'Haplotypes', table([
        row('Alignment', f'QC-passed sequences aligned with MAFFT {v(tools.get("mafft"), "v{}", "")} (<code>--auto</code>)'),
        row('Core region', f'Alignment columns with sequence in at least {v(core.get("core_min_coverage"))} of '
                           f'sequences; {v(applied.get("core_region_length_bp"))} bp in this run'),
        row('Core coverage per record', f'Records whose ungapped core is shorter than {v(core.get("min_sequence_coverage"))} '
                                         f'of the 95th-percentile core length ({v(applied.get("min_sequence_length_bp"), "{:.0f}")} bp) '
                                         f'are excluded ({v(summary.get("n_excluded_core"), "{:,}")} records)'),
        row('Haplotypes', 'Each distinct core sequence is a haplotype (exact sequence variant), keyed by the md5 of '
                          'its sequence; haplotypes that differ only at gaps or Ns are merged'),
        row('Error filter', f'Singleton haplotypes within p-distance {v(d)} of another haplotype are treated as '
                            'likely sequencing errors and removed'),
        row('Correction', ('A removed singleton record is assigned to its nearest retained haplotype when that '
                           f'haplotype is unique and within {v(d)}' if hap.get('assign_corrected', True)
                           else 'Disabled (--no-assign-corrected): removed singleton records stay unassigned')),
        row('Suspect flags', f'Singletons more than {v(hap.get("max_singleton_distance"))} from their nearest neighbor, '
                             'distance outliers and members failing the reading-frame check flag a haplotype as suspect'),
        row('IDs', 'H0001, H0002, ... in order of member count, then seq_md5'),
        row('Haplotypes retained', v(applied.get('n_haplotypes'), "{:,}")),
    ]))

    status_rows = [row(s, v(summary.get(c), "{:,}")) for s, c in (
        ('assigned', 'n_assigned'), ('assigned_corrected', 'n_assigned_corrected'),
        ('excluded_error_singleton', 'n_excluded_error'), ('excluded_core_coverage', 'n_excluded_core'),
        ('failed_qc', 'n_failed_qc'), ('no_sequence', 'n_no_sequence'))]
    status_rows.append(row('Assigned (of QC-passed records)', v(summary.get('pct_assigned'), "{}%")))
    html += tab('statuses', 'Record Statuses', '<p>Every unique record has one status (definitions in README.md).</p>\n'
                + table(status_rows))

    html += tab('taxonomy', 'Taxonomy &amp; Labels', table([
        row('Record label', 'The species field, or the genus when the species is missing'),
        row('Haplotype taxon', f'The species carried by more than {v(tax.get("majority_species_threshold"))} of the '
                               'counted labels (members with a species-level label, plus genus-only or open labels '
                               'of another genus; members without a species label are left out; ties not resolved); '
                               'otherwise the most common genus'),
        row('Label conflicts', 'Members whose species-level label differs from the haplotype’s majority label '
                               '(or from its curated taxon); listed with their evidence in label_conflicts.csv'),
        row('Divergence', 'Uncorrected p-distance between aligned haplotype sequences, skipping gaps and Ns; '
                          'no taxonomic thresholds are applied'),
    ] + ([row('Curated taxa', f'{v(curation.get("n_applied"))} haplotype(s) given the taxon in the curation file '
                              f'<code>{html_escape(str(curation.get("path", "")))}</code> (md5 {v(curation.get("md5"))}); '
                              f'{v(curation.get("n_not_applied"))} row(s) not applied (run_provenance.json)')]
         if curation else [])))

    if tree_built:
        body = table([
            row('Tips', f'Haplotypes with at least {v(applied.get("tree_min_members"))} member(s) and at least '
                        f'{v(applied.get("tree_min_length_bp"))} bp, not suspect: {v(applied.get("n_tree_tips"))} tips, '
                        f'{v(applied.get("n_tree_excluded"))} excluded (phylogeny/excluded_tips.txt)'),
            row('Alignment', f'MAFFT; trimmed with trimAl ({v(phy.get("trim_method"))}) when trimAl is installed'
                if phy.get('trim_alignment') else 'MAFFT'),
            row('Inference', f'FastTree {v(tools.get("fasttree"), "v{}", "")} <code>-nt -gtr -gamma</code>, '
                             f'seed {v(prov.get("seeds", {}).get("fasttree"))}; support values are FastTree local support'),
            row('Rooting', rooting_text),
        ])
    else:
        body = '<p class="alert alert-info">No tree was built in this run (--build-tree not given).</p>\n'
    html += tab('phylogenetic', 'Phylogeny', body)

    if geo:
        ref = (f'Custom shapefile <code>{Path(custom_shp).name}</code>, field {v(applied.get("shapefile_field"))}'
               if custom_shp else 'Global Oceans and Seas (GOaS) v1')
        body = table([
            row('Regions', f'{ref}; region column <code>{geo_cat}</code>'),
            row('Records assigned a region', 'Only records whose coordinates are of geographic quality '
                                             '(coord_status ok). Country- and province-centroid coordinates are '
                                             'kept but not assigned a region or mapped'),
            row('Coordinate status counts', ', '.join(f'{k} {n:,}' for k, n in notes.get('coord_status', {}).items())),
        ])
    else:
        body = '<p class="alert alert-info">No geographic analysis in this run.</p>\n'
    html += tab('geographic', 'Geography', body)

    sw = [row('BOLDGenotyper', f'v{prov.get("boldgenotyper_version", version)}'),
          row('Python', v(prov.get('python_version')))]
    sw += [row(name, v(tools.get(key))) for name, key in (('MAFFT', 'mafft'), ('FastTree', 'fasttree'),
                                                            ('trimAl', 'trimal')) if tools.get(key)]
    html += tab('software', 'Software', table(sw))

    n_core_ok = summary.get('n_core_ok')
    text = (f"COI barcode records for {organism} were downloaded from the Barcode of Life Data System "
            f"(BOLD; Ratnasingham &amp; Hebert, 2007) and analyzed with BOLDGenotyper v{prov.get('boldgenotyper_version', version)}. "
            f"Of {v(summary.get('n_unique'), '{:,}', '[N]')} unique records, {v(summary.get('n_qc_pass'), '{:,}', '[N]')} "
            f"passed quality control (length, ambiguous bases and an open reading frame under genetic code {v(code, missing='[N]')}). "
            f"Sequences were aligned with MAFFT (Katoh &amp; Standley, 2013), and the {v(applied.get('core_region_length_bp'), missing='[N]')} bp "
            f"core region was extracted; {v(n_core_ok, '{:,}', '[N]')} records covered it. Each distinct core sequence was treated "
            f"as a haplotype; singleton haplotypes within p-distance {v(d, missing='[d]')} of another haplotype were treated as "
            f"likely sequencing errors, and their records assigned to the nearest retained haplotype when unique and within "
            f"that distance. {v(summary.get('n_assigned', 0) + summary.get('n_assigned_corrected', 0) if summary else None, '{:,}', '[N]')} "
            f"records ({v(summary.get('pct_assigned'), missing='[X]')}% of QC-passed records) were assigned to "
            f"{v(applied.get('n_haplotypes'), '{:,}', '[M]')} haplotypes. ")
    if tree_built:
        text += ("A tree of the haplotypes was inferred with FastTree (Price et al., 2010) under the GTR+Gamma model"
                 + {'midpoint': " and rooted at its midpoint. ", 'none': ". "}.get(
                     rooting.get('method'), " and rooted on an outgroup. " if rooting.get('method') else ". "))
    if geo and not custom_shp:
        text += ("Records with geographic-quality coordinates were assigned to ocean basins of the Global Oceans and "
                 "Seas dataset (Flanders Marine Institute, 2021). ")
    elif geo:
        text += f"Records with geographic-quality coordinates were assigned to regions of a custom shapefile ({geo_cat}). "
    html += tab('statement', 'Methods Statement',
                '<div class="code-block" style="background-color: #f8f9fa; padding: 15px; border-left: 4px solid #007bff; '
                'margin: 10px 0; font-family: \'Times New Roman\', serif; line-height: 1.6;">\n'
                f'<p>{text}</p>\n</div>\n'
                '<button onclick="navigator.clipboard.writeText(this.previousElementSibling.innerText)" class="btn" '
                'style="margin-top: 10px;">Copy to Clipboard</button>\n')

    refs = ['<p>Katoh, K., &amp; Standley, D. M. (2013). MAFFT multiple sequence alignment software version 7: improvements in performance and usability. <em>Molecular Biology and Evolution</em>, 30(4), 772-780.</p>\n',
            '<p>Ratnasingham, S., &amp; Hebert, P. D. (2007). BOLD: The Barcode of Life Data System. <em>Molecular Ecology Notes</em>, 7(3), 355-364.</p>\n']
    if tree_built:
        refs.append('<p>Price, M. N., Dehal, P. S., &amp; Arkin, A. P. (2010). FastTree 2 – Approximately maximum-likelihood trees for large alignments. <em>PLoS ONE</em>, 5(3), e9490.</p>\n')
        refs.append('<p>Bianchini, G., &amp; Sánchez-Baracaldo, P. (2024). TreeViewer: Flexible, modular software to visualise and manipulate phylogenetic trees. <em>Ecology and Evolution</em>, 14, e10873. <a href="https://doi.org/10.1002/ece3.10873" target="_blank">https://doi.org/10.1002/ece3.10873</a></p>\n')
        if tools.get('trimal') and phy.get('trim_alignment'):
            refs.append('<p>Capella-Gutiérrez, S., Silla-Martínez, J. M., &amp; Gabaldón, T. (2009). trimAl: a tool for automated alignment trimming in large-scale phylogenetic analyses. <em>Bioinformatics</em>, 25(15), 1972-1973.</p>\n')
    if geo and not custom_shp:
        refs.append('<p>Flanders Marine Institute (2021). Global Oceans and Seas, version 1. Available online at <a href="https://www.marineregions.org/" target="_blank">https://www.marineregions.org/</a></p>\n')
    html += tab('references', 'References', '<div style="font-size: 0.9em; line-height: 1.8;">\n' + ''.join(sorted(refs)) + '</div>\n')

    html += '</div>\n'
    return html


def generate_html_report(
    organism: str,
    output_dir: Path,
    version: Optional[str] = None,
    staging_dir: Optional[Path] = None,
    report_path: Optional[Path] = None,
    links: Optional[dict] = None,
    provenance: Optional[dict] = None,
) -> Optional[Path]:
    """
    Generate comprehensive HTML summary report.

    Aggregates all pipeline outputs into a single, formatted HTML document
    with tables, metrics, and analysis summaries.

    Parameters
    ----------
    organism : str
        Organism name
    output_dir : Path
        Run folder (run_provenance.json, the log, summary.csv)
    version : str, optional
        BOLDGenotyper version (default: the installed package version)
    staging_dir : Path, optional
        Folder holding the step outputs (taxonomy/, quality_control/, ...);
        default output_dir (a v1-style folder)
    report_path : Path, optional
        Where to write the report (default: output_dir/report.html)
    links : dict, optional
        Staged file -> path in the run folder (layout.plan_moves). Figures are
        linked at those paths; staged files not in it are not shown. Without
        it, figures are linked relative to the report.
    provenance : dict, optional
        The pipeline's in-memory provenance (RunProvenance.to_dict()); the
        Methods and Parameters sections read it. Read from run_provenance.json
        when omitted.

    Returns
    -------
    Optional[Path]
        Path to generated HTML report, or None if generation failed

    Examples
    --------
    >>> from pathlib import Path
    >>> report_path = generate_html_report(
    ...     organism="Sphyrna_lewini",
    ...     output_dir=Path("data/Sphyrna_lewini_output")
    ... )
    >>> print(f"Report generated: {report_path}")
    """
    logger.info(f"Generating HTML summary report for {organism}...")
    if version is None:
        from . import __version__ as version

    try:
        # Initialize builder
        builder = HTMLReportBuilder(organism=organism, version=version)

        # Define paths: step outputs come from staging, run-level files from the run folder
        run_dir = Path(output_dir)
        output_dir = Path(staging_dir) if staging_dir is not None else run_dir
        reports_dir = output_dir / 'reports'
        taxonomy_dir = output_dir / 'taxonomy'
        quality_dir = output_dir / 'quality_control'

        assignment_summary_file = run_dir / "summary.csv"
        if not assignment_summary_file.exists():
            assignment_summary_file = reports_dir / f"{organism}_assignment_summary.csv"
        annotated_file = output_dir / f"{organism}_annotated.csv"
        # Note: diagnostics may be in 'haplotype_assignments' or 'genotype_assignments' (legacy)
        diagnostics_file = output_dir / 'haplotype_assignments' / f"{organism}_diagnostics.csv"
        if not diagnostics_file.exists():
            # Try legacy path
            diagnostics_file = output_dir / 'genotype_assignments' / f"{organism}_diagnostics.csv"

        # Species composition directory (haplotype_assignments for current workflow)
        species_comp_dir = output_dir / 'haplotype_assignments'
        if not species_comp_dir.exists():
            species_comp_dir = output_dir / 'genotype_assignments'

        # Load data
        assignment_summary = pd.DataFrame()
        if assignment_summary_file.exists():
            assignment_summary = pd.read_csv(assignment_summary_file)

        annotated_df = pd.DataFrame()
        if annotated_file.exists():
            try:
                annotated_df = pd.read_csv(annotated_file, low_memory=False, quoting=1)
            except Exception as e:
                logger.warning(f"Could not load annotated CSV: {e}", extra={"step_failure": "12 HTML report"})

        diagnostics_df = pd.DataFrame()
        if diagnostics_file.exists():
            diagnostics_df = pd.read_csv(diagnostics_file)

        # Add quick stats
        if not assignment_summary.empty:
            summary = assignment_summary.iloc[0]
            builder.add_quick_stat(
                _format_number(summary.get('n_unique', 0), 0),
                'Unique Records'
            )
            builder.add_quick_stat(
                _format_number(summary.get('n_haplotypes', 0), 0),
                'Haplotypes'
            )
            builder.add_quick_stat(
                _format_percentage(summary.get('pct_assigned', 0)),
                'Assignment Rate'
            )

        # Build sections
        builder.add_section(_build_executive_summary_section(assignment_summary, annotated_df))
        builder.add_section(_build_parameters_section(run_dir, organism, provenance))
        resolve = _link_resolver(links) if links is not None else None
        builder.add_section(_build_visualizations_section(output_dir, organism, resolve))
        builder.add_section(_build_methods_section(run_dir, organism, version, provenance))
        builder.add_section(_build_assignment_section(assignment_summary, diagnostics_df))
        builder.add_section(_build_quality_control_section(quality_dir))
        builder.add_section(_build_taxonomy_section(taxonomy_dir, organism))
        builder.add_section(_build_species_composition_section(species_comp_dir, organism, annotated_df))
        builder.add_section(_build_geographic_section(annotated_df))
        builder.add_section(_build_metadata_summary_section(output_dir, organism, resolve))

        # Render HTML
        html_content = builder.render()

        # Save report
        output_file = Path(report_path) if report_path is not None else run_dir / "report.html"
        output_file.parent.mkdir(parents=True, exist_ok=True)

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        logger.info(f"HTML report saved: {output_file}")
        return output_file

    except Exception as e:
        logger.error(f"Failed to generate HTML report: {e}", exc_info=True)
        return None
