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

import logging
from typing import Dict, List

logger = logging.getLogger(__name__)

# LogRecord attribute that names the step whose failure a warning reports.
STEP_FAILURE = "step_failure"


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
