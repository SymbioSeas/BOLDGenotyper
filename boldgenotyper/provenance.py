"""Run provenance and step warnings for a BOLDGenotyper run.

StepWarnings collects failures of optional steps and figures so that a run
reports them in one place (the end of the log and run_provenance.json) instead
of only in scattered log lines.

Functions that catch a failure themselves and carry on (for example a plotting
function whose savefig fails) mark the warning they log with
``extra={STEP_FAILURE: "<step>"}``. While a StepWarnings collector is attached,
it records those marked warnings too. Unmarked warnings (data-quality notices)
stay in the log only.
"""

import hashlib
import json
import logging
import platform
import re
import shlex
import shutil
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Dict, List

logger = logging.getLogger(__name__)


def count_data_rows(tsv_path) -> int:
    """Number of non-empty lines after the header."""
    with open(tsv_path, "r", encoding="utf-8", errors="replace") as f:
        return max(sum(1 for line in f if line.strip()) - 1, 0)

# LogRecord attribute that names the step whose failure a warning reports.
STEP_FAILURE = "step_failure"

FASTTREE_SEED = 314159  # FastTree's own default, now passed explicitly and recorded

_VERSION_COMMANDS = {
    "mafft": (["mafft", "--version"], r"v?(\d+\.\d+\S*)"),
    "fasttree": (["FastTree", "-expert"], r"FastTree(?: version)? (\d+\.\d+\S*)"),
    "trimal": (["trimal", "--version"], r"trimAl v?(\S+)"),
}


class _CaptureHandler(logging.Handler):
    """Record warnings marked with STEP_FAILURE into a StepWarnings collector."""

    def __init__(self, owner):
        super().__init__(level=logging.WARNING)
        self.owner = owner

    def emit(self, record):
        step = getattr(record, STEP_FAILURE, None)
        if step:
            self.owner._items.append({"step": str(step), "message": record.getMessage()})


class StepWarnings:
    """Collect step and figure failures during a pipeline run."""

    def __init__(self):
        self._items: List[Dict[str, str]] = []
        self._handler = None
        self._attached_to = None

    def attach(self, logger_name="boldgenotyper"):
        """Start recording marked warnings logged under logger_name."""
        self.detach()
        self._handler = _CaptureHandler(self)
        self._attached_to = logging.getLogger(logger_name)
        self._attached_to.addHandler(self._handler)
        return self

    def detach(self):
        """Stop recording marked warnings."""
        if self._handler is not None:
            self._attached_to.removeHandler(self._handler)
            self._handler = None
            self._attached_to = None

    def record(self, step, message):
        """Record a failure of `step` and log it as a warning."""
        item = {"step": str(step), "message": str(message)}
        self._items.append(item)
        logger.warning(f"[{item['step']}] {item['message']}")

    def __len__(self):
        return len(self._items)

    def as_list(self):
        return [dict(i) for i in self._items]

    def log_summary(self, log, organism):
        """Log the warnings summary; the last line states the run status."""
        if self._items:
            log.info(f"WARNINGS ({len(self._items)}):")
            for i in self._items:
                log.info(f"  - [{i['step']}] {i['message']}")
            noun = "warning" if len(self._items) == 1 else "warnings"
            log.info(f"Pipeline completed with {len(self._items)} {noun} for {organism}")
        else:
            log.info(f"Pipeline completed successfully for {organism}")


def file_md5(path) -> str:
    """md5 hex digest of a file."""
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tool_versions() -> Dict[str, object]:
    """Versions of MAFFT, FastTree and trimAl on PATH (None when not installed)."""
    out = {}
    for name, (cmd, pattern) in _VERSION_COMMANDS.items():
        exe = shutil.which(cmd[0]) or shutil.which(name)
        if exe is None:
            out[name] = None
            continue
        try:
            res = subprocess.run([exe] + cmd[1:], capture_output=True, text=True,
                                 timeout=20, stdin=subprocess.DEVNULL)
            m = re.search(pattern, (res.stdout or "") + (res.stderr or ""), flags=re.IGNORECASE)
            out[name] = m.group(1) if m else "unknown"
        except (OSError, subprocess.SubprocessError):
            out[name] = "unknown"
    return out


def git_commit_info(package_dir) -> Dict[str, object]:
    """Commit and dirty flag when the package is imported from a git checkout."""
    repo = Path(package_dir).resolve().parent
    if not (repo / ".git").exists() or shutil.which("git") is None:
        return {}
    try:
        sha = subprocess.run(["git", "--no-optional-locks", "-C", str(repo), "rev-parse", "HEAD"],
                             capture_output=True, text=True, timeout=10).stdout.strip()
        dirty = subprocess.run(["git", "--no-optional-locks", "-C", str(repo), "status",
                                "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return {}
    return {"git_commit": sha, "git_dirty": bool(dirty)} if sha else {}


def _now() -> str:
    return datetime.now().astimezone().isoformat(timespec="seconds")


class RunProvenance:
    """Everything needed to reproduce and audit one pipeline run.

    Built from the fully resolved PipelineConfig (every field), plus the
    environment, the input file, values computed during the run (`applied`),
    data-quality notes, and the step warnings.
    """

    def __init__(self, cfg, argv, input_path, output_dir, download_date=None):
        from . import __version__
        from .config import _convert_paths_to_strings
        self.version = __version__
        self.config = _convert_paths_to_strings(cfg.to_dict())
        self.argv = [str(a) for a in argv]
        self.input_path = Path(input_path)
        self.output_dir = Path(output_dir)
        self.started_at = _now()
        self.finished_at = None
        self.status = "running"
        self.applied = {}
        self.data_notes = {}
        self.warnings = StepWarnings()
        self.input = {
            "path": str(self.input_path),
            "md5": file_md5(self.input_path),
            "n_rows": count_data_rows(self.input_path),
            "download_date": download_date,
        }
        self._tools = tool_versions()
        self._git = git_commit_info(Path(__file__).parent)

    def to_dict(self):
        return {
            "schema_version": 1,
            "boldgenotyper_version": self.version,
            **self._git,
            "python_version": platform.python_version(),
            "platform": platform.platform(),
            "argv": self.argv,
            "command_line": shlex.join(self.argv),
            "input": dict(self.input),
            "output_dir": str(self.output_dir),
            "config": self.config,
            "applied": self.applied,
            "tools": self._tools,
            "seeds": {"fasttree": FASTTREE_SEED},
            "started_at": self.started_at,
            "finished_at": self.finished_at,
            "status": self.status,
            "warnings": self.warnings.as_list(),
            "data_notes": self.data_notes,
        }

    def finish(self, success=True):
        """Set the end time and final status."""
        self.finished_at = _now()
        if not success:
            self.status = "failed"
        else:
            self.status = "completed_with_warnings" if len(self.warnings) else "completed"

    def write(self, path):
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=2, default=str))
        return path


_PREVIOUS_RUN_FILES = ("run_provenance.json", "records.csv", "haplotypes.fasta")
_PREVIOUS_RUN_GLOBS = ("*_pipeline_parameters.json", "*_annotated.csv")


def find_previous_run(output_dir) -> List[str]:
    """Names of files or folders in output_dir that show an earlier run.

    The current run's own log is not a marker, and neither is an empty folder
    (a run that failed early can leave empty folders behind).
    """
    out = Path(output_dir)
    if not out.is_dir():
        return []
    found = [name for name in _PREVIOUS_RUN_FILES if (out / name).exists()]
    for pattern in _PREVIOUS_RUN_GLOBS:
        found.extend(p.name for p in out.glob(pattern))
    hap = out / "haplotypes"
    if hap.is_dir() and any(hap.iterdir()):
        found.append("haplotypes/")
    return sorted(set(found))
