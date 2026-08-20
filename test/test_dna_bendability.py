# Tests for src/dna_bendability.py


def test_build_backbone_length_matches_sequence():
    """Output coords array should have one row per input base."""


def test_build_backbone_unknown_dinucleotide_uses_neutral_fallback():
    """A sequence containing a base pair not in WEDGE_PARAMS (e.g. containing
    'N') should use the (0, 0, 34.3) fallback rather than raising a KeyError."""


def test_local_curvature_nan_at_boundaries():
    """First and last `window//2` positions should be NaN; interior values
    for a straight/idealised backbone should be near zero."""


def test_local_curvature_handles_degenerate_zero_length_vectors():
    """Repeated identical coordinates (norm < 1e-10) should be skipped
    (left as NaN) rather than causing a division by zero."""


def test_writhe_sign_convention():
    """A backbone constructed to be right-handed vs left-handed supercoiled
    (or a simple planar/straight backbone with writhe ~0) should return the
    expected sign/magnitude."""


def test_bend_plane_angle_planar_backbone_is_zero():
    """A backbone confined to the XY plane should return bend_plane_angle
    close to 0 degrees."""


def test_bend_plane_angle_perpendicular_backbone_is_90():
    """A backbone confined to a plane perpendicular to XY should return
    close to 90 degrees."""


def test_pair_sequences_matches_by_stripped_rc_suffix():
    """Accession names differing only by a trailing '_rc' should still be
    matched between left_names and right_names."""


def test_pair_sequences_drops_unmatched_and_wrong_length_entries():
    """Left/right entries with no matching accession, or with length !=
    flank_len, should be silently excluded from the combined output --
    confirm analyse_dna_curvature's paired-count log line reflects the
    drop (e.g. "Paired X/Y left and Z right flanks...") so a large drop
    is at least visible, even though _pair_sequences itself logs nothing."""


def test_analyse_dna_curvature_end_to_end_smoke(tmp_path):
    """Small synthetic paired left/right sequences should produce three PNGs
    (curvature, writhe, bendplane) under <outdir>/graphs/ and a result dict
    with n_paired, mean_writhe, mean_bend_angle, without raising."""


def test_analyse_dna_curvature_empty_combined_seqs_does_not_crash():
    """If no left/right names match (n_paired == 0), mean_writhe and
    mean_bend_angle should be NaN, not raise from np.mean([]) -- verify the
    current code's guard (`if len(combined_writhes) else float('nan')`)
    actually prevents the RuntimeWarning/crash in practice."""
