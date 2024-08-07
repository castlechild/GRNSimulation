from setuptools import setup, find_packages

setup(
    name='ochunGRN',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'numpy',
        'scipy',
        'networkx',
        'matplotlib'
    ],
    extras_require={
        'dev': ['pytest'],
    },
    package_data={
        "ochunGRN": ["ODESystems/41586_2011_BFnature10098_MOESM304_ESM.xls"],
        "ochunGRN": ["GRN/41598_2021_3625_MOESM5_ESM.xlsx"],
    },
    author='Melvin Bonamour',
    author_email='melaug91@hotmail.fr',
    description="Simulation dynamique de Benchmark de Reseau de Regulation Genique pour l'inférence",
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/castlechild/GRNSimulation',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
