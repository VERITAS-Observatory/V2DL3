import logging

from pyV2DL3.VtsDataSource import VtsDataSource, VtsDataSourceMeta


class EventDisplayDataSource(VtsDataSource, metaclass=VtsDataSourceMeta, file_type="ED"):
    def __init__(self, etv_file, ea_file):
        super(EventDisplayDataSource, self).__init__("EventDisplay")
        self.__evt_file__ = etv_file
        self.__ea_file__ = ea_file
    
    def __fill_evt__(self, *args, **kwargs):
        try:
            import yaml
            
            with open(kwargs["evt_filter"], "r") as file:
                evt_filter = yaml.load(file, Loader=yaml.FullLoader)
        except (KeyError, TypeError):
            # evt_filter option not used
            evt_filter = {}
        
        self._fill_evt(evt_filter, *args, **kwargs)
    
    def __fill_response__(self, *args, **kwargs):
        logging.info(
            (
                "Parameters used to query IRFs:"
                " az={0:.2f} deg,"
                " ze={1:.2f} deg,"
                " pedvar={2:.1f},"
                " offset={3:.2f} deg"
            ).format(
                self.__azimuth__,
                self.__zenith__,
                self.__pedvar__,
                self.__offset__)
        )
        
        self._fill_response(*args, **kwargs)
