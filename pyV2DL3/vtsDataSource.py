class vtsDataSource(object):
    def __init__(self,evt_file,ea_file):
        self.__evt_file__ = evt_file
        self.__ea_file__  = ea_file
        
        self.__evt__ = dict()
        self.__git__ = dict()
        self.__response__  = dict()

    def fill_data(self,**kwargs):            
        self.__fill_evt__(**kwargs)
        self.__fill_gti__(**kwargs)
        self.__fill_response__(**kwargs)

    def __fill_evt__(self):
        pass 

    def __fill_gti__(self):
        pass 

    def __fill_response__(self):
        pass

    def get_evt_data(self):
        return self.__evt__

    def get_gti_data(self):
        return self.__gti__
    
    def get_response_data(self):
        return self.__response__

