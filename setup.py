from setuptools import setup, find_packages

setup(
    name='pyV2DL3',
    version='0.1',
    packages=find_packages(),
    install_requires=[
        'click',
        'numpy',
        'scipy',
        # 'gammapy',
        'astropy',
        'root_numpy',
        'pkgconfig',
        'uproot4'
    ],
    entry_points='''
        [console_scripts]
        v2dl3=pyV2DL3.script.v2dl3:cli
        v2dl3_qsub=pyV2DL3.script.v2dl3_qsub:cli
        generate_index_file=pyV2DL3.script.generate_index_file:cli
    '''
)
