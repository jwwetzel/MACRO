"""Unit tests for ``macro_phot.external`` — the external-record arithmetic.

Every test states the fact about the SCIENCE it protects, not just the fact
about the code.  Synthetic inputs whose right answer is known by
construction, and at least one test per function that would fail if the
function were quietly replaced by something plausible-but-wrong.

The four things that could silently corrupt the YZ Cnc branch decision are
each defended by name:

  1. the night label is local, not UTC   -> test_night_label_*
  2. our own data came back as "external" -> test_is_own_*, test_nightly_*
  3. an outburst is not a superoutburst   -> test_find_episodes_*,
                                             test_plateau_*
  4. a bracket argument must be physical  -> test_superoutburst_excluded_*
"""

from __future__ import annotations

import math

import pytest

from macro_phot import external as ex


# ==========================================================================
# 1 - time: the off-by-one that would move every verdict onto its neighbour
# ==========================================================================

def test_jd_to_utc_night_is_the_calendar_date_of_the_jd():
    """JD 2460361.668 is 2024-02-21 04:02 UTC.  If this returns 02-20 the
    whole state timeline slides by a day."""
    assert ex.utc_night(2460361.668654) == "2024-02-21"


def test_local_night_label_walks_an_after_midnight_frame_back_one_day():
    """The manifest's `night` is the LOCAL evening date.  A frame taken at
    04:02 UTC on Feb 21 belongs to the evening of Feb 20 in Arizona."""
    assert ex.local_night_label(2460361.668654) == "2024-02-20"


def test_night_label_to_utc_date_is_the_inverse_for_this_season():
    """The strategy quotes the dense blocks in UTC ("Feb 21-24"); the
    manifest calls them 2024-02-20..23.  These must map onto each other, or
    the classifier tags the wrong nights."""
    assert ex.night_label_to_utc_date("2024-02-20") == "2024-02-21"
    assert ex.night_label_to_utc_date("2024-05-02") == "2024-05-03"


def test_jd_mjd_round_trip():
    for jd in (2460361.5, 2400000.5, 2451545.0):
        assert ex.mjd_to_jd(ex.jd_to_mjd(jd)) == pytest.approx(jd)


def test_days_between_is_signed_and_whole():
    assert ex.days_between("2024-02-20", "2024-02-25") == 5
    assert ex.days_between("2024-02-25", "2024-02-20") == -5


# ==========================================================================
# 2 - independence: our own photometry must never count as external
# ==========================================================================

def test_is_own_observation_fires_on_the_observer_code():
    assert ex.is_own_observation("MALW", "")


def test_is_own_observation_fires_on_the_telescope_comment():
    """A different submitter code with the RLMT comment is still our data."""
    assert ex.is_own_observation(
        "XXXX", "TAKEN WITH MACRO CONSORTIUM'S ROBERT L. MUTEL TELESCOPE "
                "AT WINER OBSERVATORY")


def test_is_own_observation_is_case_insensitive():
    assert ex.is_own_observation("xxx", "taken with the mutel telescope")


def test_a_genuine_outside_observer_is_not_flagged_as_ours():
    """The test that keeps the flag from swallowing the whole record and
    leaving the verdict with no independent evidence at all."""
    assert not ex.is_own_observation("FJQ", "CCD V-band, transformed")


def test_nightly_independent_only_drops_our_rows():
    """Two observers on one night, one of them us: the independent-only
    aggregate must see exactly one point."""
    obs = [
        ex.Obs(jd=2460361.7, mag=14.0, err=None, band="V", source="aavso",
               observer="FJQ", is_own=False),
        ex.Obs(jd=2460361.8, mag=12.0, err=None, band="V", source="aavso",
               observer="MALW", is_own=True),
    ]
    both = ex.nightly(obs, bands=ex.V_LIKE_BANDS)
    ind = ex.nightly(obs, bands=ex.V_LIKE_BANDS, independent_only=True)
    assert both[0].n == 2 and both[0].n_independent == 1
    assert ind[0].n == 1
    # And the magnitude must be the independent one, not the blend.
    assert ind[0].mag == pytest.approx(14.0)


def test_nightly_uses_the_median_not_the_mean():
    """A dwarf nova flickers within a night; one bright outlier must not
    drag the night's state tag upward."""
    obs = [ex.Obs(jd=2460361.7 + i * 0.01, mag=m, err=None, band="V",
                  source="aavso") for i, m in enumerate([14.0, 14.1, 10.0])]
    assert ex.nightly(obs, bands=ex.V_LIKE_BANDS)[0].mag == pytest.approx(14.0)


def test_nightly_excludes_fainter_than_limits():
    """A non-detection is not a magnitude.  Averaging one in would invent a
    brightness the observer explicitly declined to report."""
    obs = [ex.Obs(jd=2460361.7, mag=14.0, err=None, band="V", source="a"),
           ex.Obs(jd=2460361.8, mag=15.5, err=None, band="V", source="a",
                  is_limit=True)]
    pts = ex.nightly(obs, bands=ex.V_LIKE_BANDS)
    assert pts[0].n == 1 and pts[0].mag == pytest.approx(14.0)


# ==========================================================================
# 3 - parsers
# ==========================================================================

AAVSO_HEADER = ("JD@@@mag@@@uncert@@@band@@@by@@@comCode@@@compStar1@@@"
                "compStar2@@@charts@@@comment@@@transformed@@@airmass@@@val@@@"
                "cmag@@@kmag@@@starName@@@obsAffil@@@mtype@@@adsRef@@@"
                "digitizer@@@credit@@@obsID@@@fainterThan@@@obsType@@@"
                "software@@@obsName@@@obsCountry")


def _aavso_row(jd, mag, band="V", by="FJQ", comment="", fainter="0"):
    cells = [str(jd), str(mag), "0.01", band, by, "", "", "", "", comment,
             "", "", "", "", "", "YZ CNC", "", "", "", "", "", "1",
             fainter, "CCD", "", "", ""]
    return "@@@".join(cells)


def test_parse_aavso_delim_reads_the_documented_columns():
    text = "\n".join([AAVSO_HEADER, _aavso_row(2460361.5, 13.5)])
    obs = ex.parse_aavso_delim(text)
    assert len(obs) == 1
    o = obs[0]
    assert o.jd == pytest.approx(2460361.5)
    assert o.mag == pytest.approx(13.5)
    assert o.band == "V" and o.observer == "FJQ" and not o.is_own


def test_parse_aavso_delim_marks_our_own_rows():
    text = "\n".join([AAVSO_HEADER,
                      _aavso_row(2460361.5, 13.5, band="SR", by="MALW",
                                 comment="TAKEN WITH ... MUTEL TELESCOPE ...")])
    assert ex.parse_aavso_delim(text)[0].is_own


def test_parse_aavso_delim_drops_unparseable_magnitudes():
    """A magnitude that is not a number is not a faint magnitude.  Coercing
    it would put an invented point into a state classification."""
    text = "\n".join([AAVSO_HEADER, _aavso_row(2460361.5, "n/a"),
                      _aavso_row(2460362.5, 13.5)])
    assert len(ex.parse_aavso_delim(text)) == 1


def test_parse_aavso_delim_carries_the_fainter_than_flag():
    text = "\n".join([AAVSO_HEADER,
                      _aavso_row(2460361.5, 15.5, fainter="1")])
    assert ex.parse_aavso_delim(text)[0].is_limit


def test_parse_aavso_delim_refuses_a_payload_with_no_magnitude_column():
    """An empty list would read as 'the star had no data'.  A payload whose
    shape we do not recognise must raise instead."""
    with pytest.raises(ValueError):
        ex.parse_aavso_delim("something@@@else\n1@@@2")


def test_parse_aavso_delim_on_empty_text_is_empty_not_an_error():
    assert ex.parse_aavso_delim("") == []


ZTF_HEADER = ("oid,expid,hjd,mjd,mag,magerr,catflags,filtercode,ra,dec,chi,"
              "sharp,filefracday,field,ccdid,qid,limitmag,magzp,magzprms,"
              "clrcoeff,clrcounc,exptime,airmass,programid")


def test_parse_ztf_csv_converts_mjd_to_jd():
    row = ("1,2,3.0,60361.5,15.2,0.03,0,zg,122.7,28.1,1,0,x,1,1,1,20,26,0.1,"
           "0,0,30,1.2,1")
    obs = ex.parse_ztf_csv("\n".join([ZTF_HEADER, row]))
    assert len(obs) == 1
    assert obs[0].jd == pytest.approx(60361.5 + 2400000.5)
    assert obs[0].band == "zg"


def test_parse_ztf_csv_keeps_flagged_points_but_marks_them():
    """'ZTF looked and got a flagged measurement' and 'ZTF did not look' are
    different statements, and the coverage accounting needs both."""
    row = ("1,2,3.0,60361.5,15.2,0.03,32768,zr,122.7,28.1,1,0,x,1,1,1,20,26,"
           "0.1,0,0,30,1.2,1")
    obs = ex.parse_ztf_csv("\n".join([ZTF_HEADER, row]))
    assert len(obs) == 1 and "32768" in obs[0].note


def test_parse_ztf_csv_header_only_is_empty():
    assert ex.parse_ztf_csv(ZTF_HEADER) == []


def test_parse_asassn_json_reads_results_and_flags_heliocentric_times():
    payload = {"count": 1, "results": [
        {"hjd": 2456251.03103, "camera": "bb", "mag": 12.196,
         "mag_err": 0.02, "flux": 51.6, "flux_err": 0.95}]}
    obs = ex.parse_asassn_json(payload)
    assert len(obs) == 1 and obs[0].band == "V"
    assert "heliocentric" in obs[0].note.lower()


def test_parse_asassn_json_tolerates_a_missing_results_key():
    assert ex.parse_asassn_json({}) == []


# ==========================================================================
# 4 - the ASAS-SN match: a wide cone must not admit the wrong star
# ==========================================================================

ASASSN_ROW = (
    "<tbody><tr>"
    "<td><a href='/variables/c7510de5-ed45-56c9-bec2-f1a93b70fad8'>"
    "ASASSN-V J081056.64+280832.7</a></td>"
    "<td>YZ Cnc</td><td>122.73599</td><td>28.14243</td><td>0.49</td>"
    "<td>13.37</td><td>3.07</td><td>95.93</td><td>UGSU</td></tr></tbody>")


def test_parse_asassn_search_extracts_the_row():
    rows = ex.parse_asassn_search(ASASSN_ROW)
    assert len(rows) == 1
    assert rows[0]["uuid"].startswith("c7510de5")
    assert rows[0]["sep_arcsec"] == pytest.approx(0.49)
    assert rows[0]["other_names"] == "YZ Cnc"


def test_pick_asassn_match_prefers_the_cross_identified_record():
    """The wide cone that is needed to find these polars can admit a
    neighbour.  A record that NAMES our target is the star; proximity alone
    is a guess."""
    near_impostor = dict(uuid="a" * 8 + "-1111-2222-3333-444444444444",
                         asassn_name="ASASSN-V J000", other_names="Some Other",
                         sep_arcsec=0.1, mean_vmag=15.0, amplitude=1.0,
                         var_type="EW")
    real = ex.parse_asassn_search(ASASSN_ROW)[0]
    got, why = ex.pick_asassn_match([near_impostor, real], "YZ Cnc")
    assert got is real
    assert "cross-identified" in why


def test_pick_asassn_match_refuses_a_distant_unnamed_record():
    far = dict(uuid="b" * 8 + "-1111-2222-3333-444444444444",
               asassn_name="ASASSN-V J999", other_names="", sep_arcsec=400.0,
               mean_vmag=15.0, amplitude=1.0, var_type="EW")
    got, why = ex.pick_asassn_match([far], "YZ Cnc")
    assert got is None and "not accepted" in why


def test_pick_asassn_match_on_no_candidates_is_a_clean_none():
    got, why = ex.pick_asassn_match([], "YZ Cnc")
    assert got is None and "no variable-star record" in why


# ==========================================================================
# 5 - the ladder
# ==========================================================================

def _pts(pairs):
    """Nightly points from (night, mag) pairs."""
    return [ex.NightPoint(night=n, mag=m, n=1, n_independent=1, bands=("V",),
                          observers=("X",), spread=0.0, independent=True)
            for n, m in pairs]


def test_quiescent_baseline_is_the_faint_decile_not_the_median():
    """A dwarf nova spends time in outburst; the plain median sits in the
    outburst tail and would understate quiescence."""
    pts = _pts([(f"2024-01-{d:02d}", m) for d, m in
                zip(range(1, 21), [15.0] * 10 + [11.0] * 10)])
    assert ex.quiescent_baseline(pts) == pytest.approx(15.0)


def test_quiescent_baseline_refuses_to_speak_from_too_few_nights():
    """A baseline from two points is a guess wearing a statistic's clothes."""
    assert ex.quiescent_baseline(_pts([("2024-01-01", 15.0),
                                       ("2024-01-02", 15.1)])) is None


def test_classify_night_has_three_outcomes_not_two():
    """A star caught mid-rise is genuinely neither quiescent nor in
    outburst, and rounding it to a neighbour would fabricate certainty."""
    base = 15.0
    assert ex.classify_night(14.9, base) == ex.STATE_QUIESCENT
    assert ex.classify_night(14.2, base) == ex.STATE_ELEVATED
    assert ex.classify_night(13.0, base) == ex.STATE_OUTBURST


def test_classify_night_on_a_missing_magnitude_says_so():
    assert ex.classify_night(None, 15.0) == ex.STATE_UNKNOWN
    assert ex.classify_night(float("nan"), 15.0) == ex.STATE_UNKNOWN


def test_band_offset_is_measured_on_shared_nights_only():
    """Placing a Sloan magnitude on a V ladder must be a measurement.  Here
    SR is exactly 0.20 fainter than V on every shared night."""
    obs = []
    for i, night_jd in enumerate([2460361.7, 2460362.7, 2460363.7]):
        obs.append(ex.Obs(jd=night_jd, mag=14.0 + i, err=None, band="V",
                          source="a"))
        obs.append(ex.Obs(jd=night_jd + 0.01, mag=14.2 + i, err=None,
                          band="SR", source="a"))
    off, n, scatter = ex.band_offset(obs, "SR")
    assert off == pytest.approx(0.20, abs=1e-9)
    assert n == 3 and scatter == pytest.approx(0.0, abs=1e-9)


def test_band_offset_returns_none_when_no_night_carries_both():
    """The band cannot be placed on the ladder, and 'None' is the answer —
    not a zero offset, which would silently assume SR == V."""
    obs = [ex.Obs(jd=2460361.7, mag=14.0, err=None, band="V", source="a"),
           ex.Obs(jd=2460371.7, mag=14.2, err=None, band="SR", source="a")]
    off, n, _ = ex.band_offset(obs, "SR")
    assert off is None and n == 0


# ==========================================================================
# 6 - episodes: an outburst is not a superoutburst
# ==========================================================================

def test_plateau_span_is_bounded_by_the_episode_it_belongs_to():
    """REGRESSION.  Walking the whole season from a faint peak reported an
    887-day 'plateau' for a one-night outburst, because nearly every night
    of a dwarf nova's record is within 1.5 mag of V=13.5.  A plateau is part
    of an event, so it is measured on the event's own nights."""
    run = _pts([("2024-03-28", 11.0), ("2024-03-29", 11.1),
                ("2024-03-30", 11.2)])
    start, end, days = ex.plateau_span(run, "2024-03-28")
    assert (start, end) == ("2024-03-28", "2024-03-30")
    assert days == pytest.approx(3.0)


def test_find_episodes_does_not_let_the_plateau_escape_into_quiet_nights():
    """REGRESSION, at the level that actually broke.  A one-night outburst
    at V=13.5 sits in a season of quiescent nights at V=14.6 — all of which
    are within PLATEAU_DEPTH_MAG of it.  Measuring the plateau against the
    SEASON walks straight through them and reports months; measuring it
    against the EPISODE gives the one night that happened."""
    base = 15.0
    season = ([(f"2024-01-{d:02d}", 14.6) for d in range(1, 15)]
              + [("2024-01-15", 13.5)]
              + [(f"2024-02-{d:02d}", 14.6) for d in range(1, 15)])
    eps = ex.find_episodes(_pts(season), base)
    assert len(eps) == 1
    assert eps[0].plateau_d == pytest.approx(1.0)
    assert eps[0].kind == ex.EPISODE_NORMAL


def test_plateau_span_stops_at_the_first_night_that_has_faded():
    run = _pts([("2024-03-28", 11.0), ("2024-03-29", 11.2),
                ("2024-03-30", 13.5), ("2024-03-31", 11.3)])
    start, end, days = ex.plateau_span(run, "2024-03-28")
    # 03-30 is 2.5 mag down: the walk stops there and never reaches 03-31,
    # which is what keeps a rebrightening out of the plateau.
    assert end == "2024-03-29" and days == pytest.approx(2.0)


def test_find_episodes_grades_a_bright_long_plateau_as_a_superoutburst():
    base = 15.0
    pts = _pts([(f"2024-03-{d:02d}", 11.0) for d in range(20, 32)])
    eps = ex.find_episodes(pts, base)
    assert len(eps) == 1
    assert eps[0].kind == ex.EPISODE_SUPEROUTBURST
    assert eps[0].peak_amp == pytest.approx(4.0)


def test_find_episodes_will_not_promote_a_short_bright_outburst():
    """Amplitude alone must not make a superoutburst: a well-caught normal
    outburst peak can be bright."""
    base = 15.0
    pts = _pts([("2024-02-22", 11.0), ("2024-02-23", 11.1),
                ("2024-02-24", 11.2)])
    eps = ex.find_episodes(pts, base)
    assert eps[0].kind == ex.EPISODE_NORMAL
    assert "plateau only" in eps[0].why


def test_find_episodes_will_not_promote_a_long_faint_excursion():
    """Duration alone must not make a superoutburst either."""
    base = 15.0
    pts = _pts([(f"2024-03-{d:02d}", 13.5) for d in range(10, 26)])
    eps = ex.find_episodes(pts, base)
    assert eps[0].kind == ex.EPISODE_NORMAL
    assert "peak only" in eps[0].why


def test_find_episodes_grades_on_the_plateau_not_the_whole_excursion():
    """REGRESSION.  A superoutburst followed by a long slow decline that
    stays above the outburst threshold produced a '38-day superoutburst'.
    The excursion may be 38 days; the PLATEAU is what grades it, and both
    are reported."""
    base = 15.0
    # 5 nights near peak, then a long tail that stays 1.2 mag above base.
    peak = [(f"2024-03-{d:02d}", 11.0) for d in range(25, 30)]
    tail = [(f"2024-04-{d:02d}", 13.7) for d in range(1, 26)]
    eps = ex.find_episodes(_pts(peak + tail), base)
    assert len(eps) == 1
    assert eps[0].duration_d > 30          # the excursion really is long
    assert eps[0].plateau_d == pytest.approx(5.0)
    assert eps[0].kind == ex.EPISODE_NORMAL


def test_find_episodes_splits_on_a_long_gap():
    """Three days of silence is long enough that joining across it would be
    invention rather than interpolation."""
    base = 15.0
    pts = _pts([("2024-01-01", 12.0), ("2024-01-02", 12.0),
                ("2024-02-01", 12.0)])
    assert len(ex.find_episodes(pts, base)) == 2


def test_find_episodes_on_a_quiet_season_is_empty():
    assert ex.find_episodes(_pts([("2024-01-01", 15.0)]), 15.0) == []


# ==========================================================================
# 7 - the bracket argument, which is what makes the verdict independent
# ==========================================================================

def test_superoutburst_excluded_when_close_faint_brackets_straddle_the_gap():
    """THE argument that keeps the Feb 21-24 verdict independent of our own
    data: a >= 8 d plateau cannot fit inside a 5 d gap whose endpoints are
    both near quiescence."""
    before = _pts([("2024-02-20", 13.85)])[0]
    after = _pts([("2024-02-25", 13.63)])[0]
    ok, why = ex.superoutburst_excluded_by_brackets(before, after, 14.65)
    assert ok and "cannot fit" in why


def test_superoutburst_not_excluded_when_the_gap_is_wide_enough_to_hide_one():
    before = _pts([("2024-02-01", 13.85)])[0]
    after = _pts([("2024-03-01", 13.63)])[0]
    ok, why = ex.superoutburst_excluded_by_brackets(before, after, 14.65)
    assert not ok and "could hide" in why


def test_superoutburst_not_excluded_when_a_bracket_is_itself_in_outburst():
    """If an endpoint is at superoutburst level the gap is INSIDE an event,
    not outside one — the exact case where the argument must not fire."""
    before = _pts([("2024-03-28", 11.0)])[0]
    after = _pts([("2024-03-31", 11.2)])[0]
    ok, why = ex.superoutburst_excluded_by_brackets(before, after, 14.65)
    assert not ok and "inside an event" in why


def test_superoutburst_not_excluded_with_a_missing_bracket():
    after = _pts([("2024-02-25", 13.63)])[0]
    ok, _ = ex.superoutburst_excluded_by_brackets(None, after, 14.65)
    assert not ok


# ==========================================================================
# 8 - the branch rule itself
# ==========================================================================

def _verdict(night, frames, state, episode="", amp=None):
    return ex.NightVerdict(utc_night=night, local_night=night,
                           n_frames=frames, state=state, mag=None, amp=amp,
                           basis="independent", evidence="", episode=episode)


def test_branch_is_superhump_only_when_a_DENSE_run_sits_in_a_superoutburst():
    v = [_verdict("2024-03-29", 300, ex.STATE_OUTBURST,
                  ex.EPISODE_SUPEROUTBURST, amp=3.6)]
    branch, why = ex.branch_recommendation(v, 90)
    assert branch == "SUPERHUMP" and "2024-03-29" in why


def test_a_superoutburst_caught_only_by_a_SHORT_run_does_not_open_the_branch():
    """Three frames inside a superoutburst is not a period analysis, and the
    rule must not be satisfiable by a snapshot."""
    v = [_verdict("2024-04-17", 3, ex.STATE_OUTBURST,
                  ex.EPISODE_SUPEROUTBURST, amp=3.0),
         _verdict("2024-02-23", 130, ex.STATE_OUTBURST, "", amp=1.86)]
    branch, _ = ex.branch_recommendation(v, 90)
    assert branch == "FALLBACK"


def test_branch_counts_dense_runs_by_state_not_by_episode_membership():
    """REGRESSION.  Counting by episode reported four OUTBURST nights as
    'elevated' purely because no outside observer covered them, so the
    summary contradicted the per-night table beneath it."""
    v = [_verdict("2024-02-22", 241, ex.STATE_OUTBURST, "", amp=1.6),
         _verdict("2024-02-23", 130, ex.STATE_OUTBURST, "", amp=1.86),
         _verdict("2024-05-03", 237, ex.STATE_QUIESCENT, "", amp=0.1)]
    branch, why = ex.branch_recommendation(v, 90)
    assert branch == "FALLBACK"
    assert "2 in outburst" in why and "0 elevated" in why


def test_branch_is_undecided_when_no_run_is_dense_enough():
    v = [_verdict("2024-03-26", 2, ex.STATE_UNKNOWN)]
    branch, _ = ex.branch_recommendation(v, 90)
    assert branch == "UNDECIDED"


# ==========================================================================
# 9 - coverage accounting and survey limits
# ==========================================================================

def test_coverage_reports_span_bands_and_a_median_night_gap():
    obs = [ex.Obs(jd=2460361.5 + 2 * i, mag=14.0, err=None,
                  band="zg" if i % 2 else "zr", source="ztf")
           for i in range(5)]
    cov = ex.coverage("yzcnc", "ztf", obs)
    assert cov.n_points == 5 and cov.n_nights == 5
    assert cov.bands == ("zg", "zr")
    assert cov.median_gap_d == pytest.approx(2.0)
    assert cov.span_d == pytest.approx(8.0)


def test_coverage_of_nothing_says_so_rather_than_inventing_a_span():
    cov = ex.coverage("stlmi", "asassn", [])
    assert cov.n_points == 0 and cov.span_d is None
    assert "no rows" in cov.notes


def test_saturation_risk_fires_for_yz_cnc_in_outburst_on_ztf():
    """YZ Cnc reaches V=10.5.  If this returned False, ZTF's silence during
    outbursts would be read as 'the star was not in outburst'."""
    at_risk, why = ex.saturation_risk("ztf", 10.5)
    assert at_risk and "bright limit" in why


def test_saturation_risk_is_quiet_for_a_faint_polar():
    at_risk, _ = ex.saturation_risk("ztf", 16.45)
    assert not at_risk


def test_aavso_has_no_bright_limit_because_observers_change_equipment():
    at_risk, _ = ex.saturation_risk("aavso", 10.5)
    assert not at_risk
