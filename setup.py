from setuptools import setup, find_packages

setup(
    name='schicluster',
    version='0.1', 

    author='Jingtian Zhou (Modified by Quanyi Zhao)',
    author_email='jiz509@eng.ucsd.edu',
    
    packages=find_packages(),
    
    description='Modified version of schicluster for migration.',
    license='MIT',
    url='https://github.com/zhoujt1994/scHiCluster',
    
    include_package_data=True,
    
install_requires=[
    'numpy>=2.0.1',
    'scipy>=1.15.3',
    'scikit-learn>=1.8.0',
    'h5py>=3.15.1',
    'joblib>=1.5.3',
    'clodius>=0.20.4',
    'tables>=3.10.2',
    'cooler>=0.10.4',
    'pandas>=2.3.3',
    'statsmodels>=0.14.6',
    'rpy2>=3.6.4',
    'anndata>=0.12.7',
    'xarray>=2025.12.0',
    'zarr>=2.18.7',
    'numcodecs>=0.15.1',
],
    
    package_data={
        '': ['*.txt', '*.tsv', '*.csv', '*.fa', '*Snakefile', '*ipynb', '*R']
    },
    
    entry_points={
        'console_scripts': [
            'hicluster=schicluster.__main__:main',
            'hic-internal=schicluster._hicluster_internal:internal_main'
        ],
    }
)