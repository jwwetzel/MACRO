"""Shared pytest bootstrap for the pipeline test suite.

The packages under test (``macro_core``, ``macro_phot``, ``rlmt_diagnostics``,
and the ``scripts`` modules) live under ``pipeline/``, which is NOT a Python
package root that pytest knows about: the documented invocation runs from
the REPO root (``python -m pytest pipeline/tests -q``), where ``pipeline/``
is not on ``sys.path`` and every ``import rlmt_diagnostics`` dies at
collection — killing the WHOLE suite, not just the new file (pytest aborts
on collection errors).

Historically each test file carried its own two-line ``sys.path.insert``
shim; a file that forgot it (the S2 review caught ``test_diagnostics.py``)
broke the repo-root invocation for everyone.  This conftest is the single
fix: pytest imports ``conftest.py`` BEFORE collecting any test module in
this directory, so putting ``pipeline/`` on ``sys.path`` here makes every
current and future test file importable from any working directory.  The
per-file shims remain harmless duplicates.
"""

import sys
from pathlib import Path

# pipeline/tests/conftest.py -> parent = pipeline/tests -> parent = pipeline/
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# ---------------------------------------------------------------------------
# "Product absent" is a FAILURE when asked to be  (finding F-9 / DS.F4)
# ---------------------------------------------------------------------------
# About 130 tests in this suite check a product database or an emitted
# manuscript file, and skip when it is not on disk.  That is the right
# default for a developer on a clean checkout — and exactly wrong for a
# release check, where "0 failed" over a run that examined none of the
# evidence is a green light nobody earned.
#
# With the environment variable MACRO_REQUIRE_PRODUCTS set, a skip whose
# reason says the product is absent is turned into a failure:
#
#     MACRO_REQUIRE_PRODUCTS=1 python -m pytest pipeline/tests -q
#
# Which skips count is decided by ONE pure function,
# ``macro_core.provenance.is_absent_product_skip`` (unit-tested in
# test_provenance.py).  Skips that mean "this assertion does not apply to
# this build" or "an optional tool is missing" stay skips in every mode.
import os                                                    # noqa: E402

import pytest                                                # noqa: E402

from macro_core import provenance as _pv                     # noqa: E402


def _products_required() -> bool:
    """True when the release mode is switched on in the environment."""
    return os.environ.get(_pv.REQUIRE_PRODUCTS_ENV, "").strip().lower() \
        not in ("", "0", "false", "no")


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Rewrite an absent-product skip as a failure, in release mode only."""
    outcome = yield
    report = outcome.get_result()
    if not (_products_required() and report.skipped):
        return
    if hasattr(report, "wasxfail"):          # an xfail is not a skip
        return
    # A skip's longrepr is (path, lineno, "Skipped: <reason>").
    longrepr = report.longrepr
    reason = longrepr[2] if isinstance(longrepr, tuple) else str(longrepr)
    if _pv.is_absent_product_skip(reason):
        report.outcome = "failed"
        report.longrepr = (
            f"PRODUCT ABSENT ({_pv.REQUIRE_PRODUCTS_ENV} is set, so this is "
            f"a failure, not a skip): {reason}")
