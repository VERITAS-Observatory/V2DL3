import ROOT

from pyV2DL3.vegas.fillEVENTS_not_safe import __fillEVENTS_not_safe__
from pyV2DL3.vegas.fillRESPONSE_not_safe import __fillRESPONSE_not_safe__
from pyV2DL3.vegas.load_vegas import VEGASStatus
from pyV2DL3.VtsDataSource import VtsDataSource


class VegasDataSource(VtsDataSource):
    def __init__(
        self,
        evt_file,
        ea_files,
        bypass_fov_cut=False,
        event_class_mode=False,
        psf_king_params=None,
        reco_type=1,
        save_msw_msl=False,
        corr_EB_params=False,
        st6_configs=None
    ):
        super(VegasDataSource, self).__init__("VEGAS", evt_file, None)

        # Loading VEGAS if not already done so
        self.vegas_status = VEGASStatus()
        self.vegas_status.loadVEGAS()
        self.__evt_file__ = ROOT.VARootIO(evt_file, True)
        if not isinstance(ea_files, list):
            ea_files = [ea_files]
        self.__ea_files__ = ea_files
        self.__event_class_mode__ = event_class_mode
        self.__fov_cut__ = not bypass_fov_cut
        self.__psf_king_params__ = psf_king_params
        self.__reco_type__ = reco_type
        self.__save_msw_msl__ = save_msw_msl
        self.__corr_EB_params__ = corr_EB_params
        self.__st6_configs__ = st6_configs

        # Auxiliary storage
        self.__azimuth__ = 0
        self.__zenith__ = 0
        self.__noise__ = 0

    def __del__(self):
        """Close the root files

        These typechecks will prevent the user from having their true exception
        buried by a CPyCppyy exception on program exit.
        """
        cpy_nonestring = "<class 'CPyCppyy_NoneType'>"

        if str(type(self.__evt_file__)) != cpy_nonestring and not isinstance(
            self.__evt_file__, str
        ):
            self.__evt_file__.closeTheRootFile()

    def __fill_evt__(self):
        gti, ea_config, evt_dicts = __fillEVENTS_not_safe__(
            self.__evt_file__,
            self.__ea_files__,
            self.__irf_to_store__,
            event_class_mode=self.__event_class_mode__,
            fov_cut=self.__fov_cut__,
            reco_type=self.__reco_type__,
            save_msw_msl=self.__save_msw_msl__,
            corr_EB=self.__corr_EB_params__,
            psf_king_params=self.__psf_king_params__,
            st6_configs=self.__st6_configs__
        )
        self.__gti__ = gti
        # This is an array of dicts for each event class (array of one when not using event class mode)
        self.__evt__ = evt_dicts
        self.__azimuth__ = ea_config["azimuth"]
        self.__zenith__ = ea_config["zenith"]
        self.__noise__ = ea_config["noise"]

    def __fill_gti__(self):
        pass

    def __fill_response__(self):
        az = self.__azimuth__
        ze = self.__zenith__
        nn = self.__noise__
        response_dicts = []
        # Fill response for each event class
        for ec in self.__ea_files__:
            response_dicts.append(
                __fillRESPONSE_not_safe__(
                    ec,
                    az,
                    ze,
                    nn,
                    self.__irf_to_store__,
                    psf_king_params=self.__psf_king_params__,
                    st6_configs=self.__st6_configs__
                )
            )

        self.__response__ = response_dicts
