# DRAFT — not sent. For James to edit and send from his own account.

**To:** Winer Observatory staff / RLMT operations *(addresses to be filled in by James)*
**Cc:** MACRO observing committee
**From:** James Wetzel (Coe College)
**Subject:** RLMT re-opening: revised request (rev. 3) — please disregard the August version
**Attachments:** `2026-10_observatory_request_rev3.md` (or PDF); `ops/eruption_block/` (README + three schedule files + CSV)

---

Hello,

In August I sent a consolidated observing and calibration request for the October
re-opening. Please set it aside. When we reviewed it last week against the frames actually
in the archive, it turned out to ask for readout modes of cameras that are no longer on
the telescope, and for T CrB runs at a time of year when the star is up for about an
hour. The attached revision 3 replaces it completely. I apologise for the churn; this one
is written for the QHY600 under pyscope as the June frames show it, and its first section
lists what we believe the configuration to be so that you can correct us.

The short version of what would help most, in order:

1. **Before anyone adjusts anything on the telescope: flats in every filter (hrg and lrg
   included), biases and darks, with the camera exactly as it is.** We hold more than
   29,000 QHY frames and not one flat for them, and once the optical train is touched the
   as-found flat is gone for good. (§3 step a, §4.)
2. **Then a through-focus run** before science — the late-June frames look astigmatic
   and about a thousand focuser counts from focus on some nights. (§3 step b.)
3. **Flat pairs at six light levels per filter.** One set of frames gives us the flat,
   the gain, the read noise and the linearity curve. (§4 item 3.)
4. **A handful of pyscope header changes** — readout mode, IMAGETYP, the gain cards, the
   filter-wheel list, sub-second DATE-OBS. (§5.)
5. **A fifteen-minute T CrB block each clear evening** until the star is lost in late
   October, with two nearby standards, resuming in the morning sky from mid-December.
   No long runs before January. (§6.A.)
6. **Eight questions** (§8). The two that matter most:
   - Do the ZWO ASI camera (2024-12 to 2026-03) and the SBIG AC4040 (2023 to 2024-03)
     still exist? A few hours of bench darks with the ASI would calibrate every T CrB
     spectrum we took in 2025.
   - Could we have a copy of the calibration directories on the control computers,
     masters and raw frames, as they stand?
7. **The T CrB eruption block:** a draft is attached. Could you load it, run it once on
   any bright star, and send the log? (§7.)

Everything else — Be-star monitoring, ST LMi from November, standards — is queue work at
whatever priority suits the consortium, and is laid out with the dates the sky allows.

If any of this conflicts with work you already have planned for the re-opening
(collimation, a camera swap, a pyscope upgrade), please tell me the order you intend and
I will fit the request around it; the only thing I would ask you to hold to is item 1
coming first. And if it is easier to talk it through, I am glad to call.

Thank you for all of it.

James Wetzel
Coe College, for the MACRO archival-analysis programme

---

*Notes for James (delete before sending):*
- *Committee source: SYNTHESIS U10 and §6 ("Send ops request rev. 3 to Winer; ask whether
  the ASI and AC4040 cameras still exist; ask for the server `calibrations/` trees and a
  hardware change log").*
- *The request file is marked DRAFT in its header; change the status line when you send.*
- *The numbers in this note (29,000 frames, no flats; fifteen minutes; late October /
  mid-December) are quoted from the generated tables in rev. 3 §2 and §4 — re-check them
  against the file if it has been regenerated since 2026-10-03.*
- *Re-opening may already have happened by the time this goes out: if science frames have
  been taken since, say so and ask for item 1 anyway — a flat taken now is still worth
  more than none.*
- *Consider converting the Markdown to PDF so the tables and the two figures travel.*
