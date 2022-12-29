from setuptools import find_packages, setup

__version__ = "1.4"

# Load README
with open('README.md', encoding='utf-8') as f:
    long_description = f.read()

setup(
    name='tclint',
    version='1.4',
    author='Josef Kehrein',
    author_email='josef.kehrein@helsinki.fi',
    license='MIT',
    description='Predict interaction of molecule with lecithin/taurocholate',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/juppifluppi/tclint',
    download_url=f'https://github.com/juppifluppi/tclint/archive/refs/tags/v_{__version__}.tar.gz',
    project_urls={
        'Source': 'https://github.com/juppifluppi/tclint',
        'Web application': 'http://tclint.streamlit.app',
    },
    packages=find_packages(),
    package_data={'tclint': ['py.typed'],
    '': ['csv.typed']},
    include_package_data=True,
    entry_points={
        'console_scripts': [
            'tclint=tclint',
        ]
    },
    install_requires=[
        'numpy>=1.23.4',        
        'dimorphite-dl>=1.3.2',
        'rdkit>=2022.09.1',    
        'scopy>=1.2.5'         
    ],
    python_requires='>=3.7',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent'
    ],
    keywords=[
        'chemistry',
        'property prediction',
        'logistic regression'
    ]
)

