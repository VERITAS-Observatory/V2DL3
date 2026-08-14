import logging


class VtsDataSource(object):
    def __init__(self, source_name, evt_file, ea_file):
        self.__data_source_name__ = source_name
        logging.info("Reconstruction type: %s", source_name)
        self.__evt_file__ = evt_file
        self.__ea_file__ = ea_file

        self.__evt__ = {}
        self.__gti__ = {}
        self.__response__ = dict()
        # 'point-like' or 'full-enclosure'
        # Default is point like
        self.__irf_to_store__ = {"point-like": True, "full-enclosure": False}

    def fill_data(self, **kwargs):
        self.__fill_evt__(**kwargs)
        self.__fill_gti__(**kwargs)
        self.__fill_response__(**kwargs)

    def __fill_evt__(self):
        pass

    def __fill_gti__(self):
        pass

    def __fill_response__(self):
        pass

    def get_source_name(self):
        return self.__data_source_name__

    def get_evt_data(self):
        return self.__evt__

    def get_gti_data(self):
        return self.__gti__

    def get_response_data(self):
        return self.__response__

    def get_version(self):
        return "0.0.0"

    def set_irfs_to_store(self, irf_to_store):
        self.__irf_to_store__ = irf_to_store
