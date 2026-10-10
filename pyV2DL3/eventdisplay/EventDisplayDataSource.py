import logging
import re

import uproot
import yaml

from pyV2DL3.eventdisplay.fillEVENTS import __fillEVENTS__
from pyV2DL3.eventdisplay.fillRESPONSE import __fill_response__
from pyV2DL3.eventdisplay.util import get_root_log_lines
from pyV2DL3.VtsDataSource import VtsDataSource


class EventDisplayDataSource(VtsDataSource):
    """
    Eventdisplay data source class holding events and IRFs.

    """
    def __init__(self, etv_file, ea_file):
        super(EventDisplayDataSource, self).__init__("EventDisplay", etv_file, ea_file)
        self.__evt_file__ = etv_file
        self.__ea_file__ = ea_file

        # Auxiliary storage
        self.__azimuth__ = 0
        self.__gti__ = None
        self.__zenith__ = 0
        self.__pedvar__ = 0

    def get_version(self):
        return self._fill_data_source_version()

    def __fill_evt__(self, **kwargs):
        evt_filter = self.__fill_event_filter(kwargs.get("evt_filter", None))

        gti, ea_config, events = __fillEVENTS__(
            self.__evt_file__, evt_filter, kwargs.get("db_fits_file", None))
        self.__gti__ = gti
        self.__evt__ = events
        self.__azimuth__ = ea_config["azimuth"]
        self.__zenith__ = ea_config["zenith"]
        self.__pedvar__ = ea_config["pedvar"]

    @staticmethod
    def __fill_event_filter(filter_file):
        """
        Read event filter from yaml file

        """
        try:
            with open(filter_file, "r", encoding="utf-8") as file:
                evt_filter = yaml.load(file, Loader=yaml.FullLoader)
        except (KeyError, TypeError):
            evt_filter = {}
        return evt_filter

    def __fill_gti__(self, **kwargs):
        pass

    def __fill_response__(self, **kwargs):
        logging.info(
            "Parameters used to query IRFs: az=%0.2f deg,"
            " ze=%1.2f deg,"
            " pedvar=%2.1f",
            self.__azimuth__,
            self.__zenith__,
            self.__pedvar__,
        )

        self.__response__ = __fill_response__(
            self.__evt_file__,
            self.__ea_file__,
            self.__azimuth__,
            self.__zenith__,
            self.__pedvar__,
            self.__irf_to_store__,
            **kwargs
        )

    def _fill_data_source_version(self):
        """
        Get Eventdisplay version from log entry in anasum file.

        Returns
        -------
        str
            Version string of the Eventdisplay data source, e.g. "v490.7"
        """
        with uproot.open(self.__evt_file__) as file:
            try:
                log = file["anasumLog;1"]
            except uproot.exceptions.KeyInFileError:
                logging.warning(f"No anasum log found in {self.__evt_file__}")
                return "0.0.0"

            try:
                data_list = get_root_log_lines(log)
            except ValueError as error:
                raise ValueError(
                    f"Could not decode anasum log in {self.__evt_file__}: {error}"
                ) from error

        matches = [line for line in data_list if "VERITAS Analysis Summary" in line]
        if len(matches) != 1:
            raise ValueError(
                f"Expected one VERITAS Analysis Summary line in "
                f"{self.__evt_file__}, found {len(matches)}"
            )
        match = re.search(r"version\s+([^\)]+)", matches[0])
        if match is None:
            raise ValueError(
                f"Could not parse Eventdisplay version from {self.__evt_file__}"
            )
        return match.group(1)
