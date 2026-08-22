from unittest.mock import patch

from pyV2DL3.script.eventdisplay_query_runparameters import get_epoch_effective_area
from pyV2DL3.tests.helpers import RootFile, RootLog


def test_query_runparameters_accepts_repeated_effective_area_lines():
    root_file = RootFile(
        {
            "anasumLog;1": RootLog(
                [
                    "reading effective areas from /aux/effArea-v490.root",
                    "reading effective areas from /aux/effArea-v490.root",
                ]
            ),
            "run_64080;1/stereo/mscwTableLog;1": RootLog(
                ["Evaluating instrument epoch (was: old, is: V6_2012_2013a)"]
            ),
        }
    )

    with patch("pyV2DL3.script.eventdisplay_query_runparameters.up.open", return_value=root_file):
        assert get_epoch_effective_area("run.anasum.root", 64080) == (
            "V6_2012_2013a",
            "effArea-v490.root",
        )
