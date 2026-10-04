def event_file(tmp_path):
    from astropy.io import fits
    from pyV2DL3.generateObsHduIndex import gen_obs_index
    names = ["OBS_ID", "RA_PNT", "DEC_PNT", "ALT_PNT", "AZ_PNT", "ONTIME", "LIVETIME", "DEADC", "TSTART", "TSTOP", "N_TELS", "NSBLEVEL", "RA_OBJ", "DEC_OBJ"]
    hdu = fits.BinTableHDU(name="EVENTS")
    for name in names: hdu.header[name] = 1
    hdu.header["TSTART"] = 683868345.678
    hdu.header["TSTOP"] = 683868346.789
    for name in ["TELLIST", "OBJECT", "DATE-OBS", "DATE-AVG", "DATE-END"]: hdu.header[name] = "value"
    for name in ["MJDREFI", "MJDREFF", "ALTITUDE", "GEOLAT", "GEOLON"]: hdu.header[name] = 0
    for name in ["TIMEUNIT", "TIMESYS", "TIMEREF"]: hdu.header[name] = "value"
    path = tmp_path / "events.fits"
    fits.HDUList([fits.PrimaryHDU(), hdu]).writeto(path)
    return path

def test_index_preserves_absolute_time_precision(tmp_path):
    from astropy.io import fits
    from pyV2DL3.generateObsHduIndex import gen_obs_index, _default_null_value
    index = gen_obs_index([event_file(tmp_path)])
    out = tmp_path / "index.fits"
    fits.HDUList([fits.PrimaryHDU(), index]).writeto(out)
    with fits.open(out) as hdul:
        assert abs(hdul[1].data["TSTART"][0] - 683868345.678) < 1e-6
        assert abs(hdul[1].data["TSTOP"][0] - 683868346.789) < 1e-6
    assert _default_null_value(">f8") == -9999.
