from astropy.io import fits
from pyV2DL3.generateObsHduIndex import gen_hdu_index

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
