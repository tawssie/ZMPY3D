from setuptools import setup, find_packages

setup(
    name='ZMPY3D',
    version='0.0.2',
    author='Jhih Siang (Sean) Lai',
    author_email='js.lai@uqconnect.edu.au, jsl035@ucsd.edu',
    description='ZMPY3D NumPy version',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/tawssie/ZMPY3D',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],

    entry_points={
        'console_scripts': [
            'ZMPY3D_CLI_ZM=ZMPY3D.ZMPY3D_CLI_ZM:main',
            'ZMPY3D_CLI_SuperA2B=ZMPY3D.ZMPY3D_CLI_SuperA2B:main',
            'ZMPY3D_CLI_ShapeScore=ZMPY3D.ZMPY3D_CLI_ShapeScore:main',
            'ZMPY3D_CLI_BatchSuperA2B=ZMPY3D.ZMPY3D_CLI_BatchSuperA2B:main',
            'ZMPY3D_CLI_BatchShapeScore=ZMPY3D.ZMPY3D_CLI_BatchShapeScore:main',
            'ZMPY3D_CLI_BatchZM=ZMPY3D.ZMPY3D_CLI_BatchZM:main',
        ],
    },

    python_requires='>=3.9.16',
    install_requires=[
        'numpy>=1.23.5',
    ],
    
    include_package_data=True, # for the cache
)



