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
