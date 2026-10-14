from astropy.io import fits

from pyV2DL3.generateObsHduIndex import gen_hdu_index, gen_obs_index


def test_long_relative_filename_and_directory_round_trip(tmp_path):
    folder = tmp_path / ("archive_" * 8)
    folder.mkdir()
    filename = "observation_" * 8 + ".fits"
    path = folder / filename
    hdu = fits.BinTableHDU(name="EVENTS")
    hdu.header.update(OBS_ID=123, HDUCLAS1="EVENTS")
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)
    index = gen_hdu_index([path], index_file_dir=tmp_path)
    out = tmp_path / "index.fits"
    fits.HDUList([fits.PrimaryHDU(), index]).writeto(out)
    with fits.open(out) as hdul:
        row = hdul[1].data[0]
        assert tmp_path / row["FILE_DIR"] / row["FILE_NAME"] == path


def test_observation_index_preserves_long_standard_and_dqm_strings(tmp_path):
    source = tmp_path / "observation.fits"
    long_object = "VERITAS source with a long catalogue designation"
    long_tellist = "T1,T2,T3,T4 and auxiliary telescope metadata"
    long_dqm_reason = "extended DQM quality explanation for this run"

    events = fits.BinTableHDU.from_columns([], name="EVENTS")
    events.header.update(
        OBS_ID=42,
        RA_PNT=83.63,
        DEC_PNT=22.01,
        ALT_PNT=70.0,
        AZ_PNT=180.0,
        ONTIME=1200.0,
        LIVETIME=1180.0,
        DEADC=0.98,
        TSTART=100.0,
        TSTOP=1300.0,
        N_TELS=4,
        TELLIST=long_tellist,
        OBJECT=long_object,
        RA_OBJ=83.63,
        DEC_OBJ=22.01,
        **{
            "DATE-OBS": "2026-10-11T00:00:00",
            "DATE-AVG": "2026-10-11T00:10:00",
            "DATE-END": "2026-10-11T00:20:00",
            "NSBLEVEL": 5.2,
            "CONFIG": 3,
            "DQMREAS": long_dqm_reason,
            "MJDREFI": 55197,
            "MJDREFF": 0.0,
            "TIMEUNIT": "s",
            "TIMESYS": "utc",
            "TIMEREF": "topocenter",
            "ALTITUDE": 1268.0,
            "GEOLAT": 31.675,
            "GEOLON": -110.952,
        },
    )
    fits.HDUList([fits.PrimaryHDU(), events]).writeto(source)

    index = gen_obs_index([str(source)], str(tmp_path), dqm_header=True)
    output = tmp_path / "obs-index.fits"
    fits.HDUList([fits.PrimaryHDU(), index]).writeto(output)

    with fits.open(output) as hdul:
        row = hdul[1].data[0]
        assert _as_text(row["OBJECT"]) == long_object
        assert _as_text(row["TELLIST"]) == long_tellist
        dqm_reason_column = next(
            name for name in hdul[1].columns.names if name.strip() == "DQMREAS"
        )
        assert _as_text(row[dqm_reason_column]) == long_dqm_reason


def _as_text(value):
    return value.decode().rstrip() if isinstance(value, bytes) else str(value).rstrip()
