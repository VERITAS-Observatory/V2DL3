from setuptools import find_packages
from setuptools import setup

setup(
    name="pyV2DL3",
    version="0.1",
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
        "root_numpy",
    ],
    entry_points="""
        [console_scripts]
        v2dl3=pyV2DL3.script.v2dl3:main
        v2dl3-vegas=pyV2DL3.script.v2dl3_for_vegas:cli
        v2dl3-eventdisplay=pyV2DL3.script.v2dl3_for_Eventdisplay:cli
        v2dl3_qsub=pyV2DL3.script.v2dl3_qsub:cli
        generate_index_file=pyV2DL3.script.generate_index_file:cli
    """,
)
