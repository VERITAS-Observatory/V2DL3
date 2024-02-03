from setuptools import find_packages
from setuptools import setup

setup(
    name="pyV2DL3",
    version="0.5",
    packages=find_packages(),
    install_requires=[
        "astropy",
        "click",
        "numpy",
        "pkgconfig",
        "pyyaml",
        "scipy",
        "tqdm",
        "uproot",
    ],
    entry_points={
        "console_scripts": [
            "v2dl3-vegas=pyV2DL3.script.v2dl3_for_vegas:cli",
            "v2dl3-generate-index-file=pyV2DL3.script.generate_index_file:cli",
            "v2dl3-eventdisplay=pyV2DL3.script.v2dl3_for_Eventdisplay:cli",
        ]
    },
)
