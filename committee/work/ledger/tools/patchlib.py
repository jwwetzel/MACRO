"""Tiny exact-anchor patch helper: every replacement must match exactly once."""
from pathlib import Path

class Patcher:
    def __init__(self, path):
        self.path = Path(path)
        self.text = self.path.read_text(encoding="utf-8")

    def rep(self, old, new, count=1):
        n = self.text.count(old)
        if n != count:
            raise SystemExit(f"anchor matched {n}x (wanted {count}):\n{old[:200]}")
        self.text = self.text.replace(old, new)

    def between(self, start, end, new):
        """Replace from `start` (inclusive) up to `end` (exclusive)."""
        i = self.text.index(start)
        if self.text.count(start) != 1:
            raise SystemExit(f"start anchor not unique: {start[:120]}")
        j = self.text.index(end, i)
        self.text = self.text[:i] + new + self.text[j:]

    def save(self):
        self.path.write_text(self.text, encoding="utf-8")
