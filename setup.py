from setuptools import setup, find_packages

setup(
    name='pyV2DL3',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'click',
        'numpy',
        'astropy',
        'root_numpy'
    ],
    entry_points='''
        [console_scripts]
        st5ToDL3=pyV2DL3.script.st5ToDL3:cli
    '''
)
