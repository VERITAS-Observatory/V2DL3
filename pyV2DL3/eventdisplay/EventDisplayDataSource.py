from pyV2DL3.VtsDataSource import VtsDataSource
from pyV2DL3.eventdisplay.fillEVENTS import __fillEVENTS__
from pyV2DL3.eventdisplay.fillRESPONSE import __fillRESPONSE__
import warnings

class EventDisplayDataSource(VtsDataSource):
    def __init__(self, etv_file, ea_file):
        super(EventDisplayDataSource, self).__init__('EventDisplay', etv_file, ea_file)
        self.__evt_file__ = etv_file
        self.__ea_file__ = ea_file

        # Auxiliary storage 
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__noise__ = 0
        self.__offset__ = 0

    def __fill_evt__(self, **kwargs):
        # can be simplified further:
        try:
            import yaml
            with open(kwargs["evt_filter"], "r") as file:
                evt_filter = yaml.load(file, Loader=yaml.FullLoader)
        except (KeyError, TypeError):
            evt_filter = {}
        except ModuleNotFoundError as e:
            warnings.warn("Failed to import yaml. Event filter will be ignored.")
            evt_filter = {}
        except FileNotFoundError:
            warnings.warn("yaml file not found. Event filter will be ignored.")
            evt_filter = {}

        gti, ea_config, events = __fillEVENTS__(self.__evt_file__, evt_filter)
        self.__gti__ = gti
        self.__evt__ = events
        self.__azimuth__ = ea_config['azimuth']
        self.__zenith__ = ea_config['zenith']
        self.__noise__ = ea_config['noise']

    def __fill_gti__(self, **kwargs):
        pass

    def __fill_response__(self, **kwargs):
        pass
        az = self.__azimuth__
        ze = self.__zenith__
        nn = self.__noise__
        oo = self.__offset__
        print("Coordinates to fillresponse:", az, ze, nn, oo)
        self.__response__ = __fillRESPONSE__(self.__evt_file__, self.__ea_file__, az, ze, nn, oo,
                                                      self.__irf_to_store__)
