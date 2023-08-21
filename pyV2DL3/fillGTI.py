from astropy.io import fits

from pyV2DL3.addHDUClassKeyword import addHDUClassKeyword
from pyV2DL3.constant import VTS_REFERENCE_MJD


def fillGTI(datasource, goodTimeStart=None, goodTimeStop=None):
    gti_dict = datasource.get_gti_data()
    goodTimeStart = gti_dict["goodTimeStart"]
    goodTimeStop = gti_dict["goodTimeStop"]
    startTime_s = gti_dict["TSTART"]
    endTime_s = gti_dict["TSTOP"]

    hdu2 = fits.BinTableHDU.from_columns(
        [
            fits.Column(name="START", format="1D", array=goodTimeStart, unit="s"),
            fits.Column(name="STOP", format="1D", array=goodTimeStop, unit="s"),
        ]
    )
    hdu2.name = "GTI"
    # Fill Standard HDUCLASS keywords
    hdu2 = addHDUClassKeyword(hdu2, class1="GTI")
    hdu2.header.set("TSTART", startTime_s, "start time [s]")
    hdu2.header.set("TSTOP", endTime_s, "stop time same [s]")

    hdu2.header.set("MJDREFI ", VTS_REFERENCE_MJD, "int part of reference MJD [days]")
    hdu2.header.set("MJDREFF ", 0.0, "fractional part of reference MJD [days]")

    hdu2.header.set("TIMEUNIT", "s", "time unit is seconds since MET start")
    hdu2.header.set("TIMESYS ", "utc", "time scale is UTC")
    hdu2.header.set("TIMEREF ", "topocenter", "location from where the observation was made")

    hdu2.header.set(
        "TIMEZERO", 0, "zero time [s]"
    )  # Times stored in GTI table are relative to reference time
    hdu2.header.set("TTYPE1", "START   ", " start of good time interval")
    hdu2.header.set("TTYPE2", "STOP    ", " start of good time interval")
    hdu2.header.set("EXTNAME", "GTI     ", " name: Good Time Intervals")

    return hdu2
