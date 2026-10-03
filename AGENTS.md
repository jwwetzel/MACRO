# MACRO Standing Committee

Seven lenses that review every plan and every finished result in this repository.
Each seat is a Claude Code subagent defined in `.claude/agents/`; spawn one by
name, or brief a general-purpose agent with "read `.claude/agents/<seat>.md` and
adopt it".

| Seat | Agent | Lens | Characteristic question |
|---|---|---|---|
| 1 | `data-scientist` | Statistics, selection effects, provenance | What is the null? |
| 2 | `observational-astronomer` | Calibration, photometry, grism practice, observing plan | Have you looked at the frames? |
| 3 | `physicist` | Physical interpretation, dimensional sanity | What does this number mean physically? |
| 4 | `detector-engineer` | Camera modes, gain, linearity, header truth | Which mode was the camera actually in? |
| 5 | `telescope-engineer` | Optics, mechanics, hardware eras, operations | What physically changed that night? |
| 6 | `journal-editor` | Venue, novelty, framing, release readiness (committee chair for manuscripts) | Why would an outside reader care? |
| 7 | `hostile-referee` | Cold adversarial read of finished work | Where, exactly, is that shown? |

## Charter

1. **Plan review.** Before a project's work begins (or resumes after a long
   gap), seats 1–6 each review the strategy, the plan ledger and the evidence on
   disk, and write a memo. The chair's synthesis becomes binding amendments to
   the project's `ANALYSIS_STRATEGY.md` and plan ledger.
2. **Result review.** When a project completes its final plan task, the result
   goes to the committee: seats 1–5 review their own domain, seat 6 reviews the
   manuscript as an editor, seat 7 writes a cold referee report.
3. **Hardening loop.** Every BLOCKER and MAJOR is fixed or rebutted in writing,
   then the seats that raised them re-review. Repeat until every seat records
   **satisfied** and seat 7's verdict is accept or minor-with-zero-majors.
4. **Independence.** A seat reviews what is on disk and is never handed the
   author's conclusion. Seat 7 reads cold.
5. **Record.** Memos live in `committee/reviews/<date>[-<project>-r<N>]/<seat>.md`;
   the chair's synthesis and the response to each finding live alongside as
   `SYNTHESIS.md` and `RESPONSE.md`. Nothing is settled orally.

Severity scale: **BLOCKER** · **MAJOR** · **MINOR** · **NOTE**.
