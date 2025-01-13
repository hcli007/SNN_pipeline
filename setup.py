from setuptools import setup, find_packages

setup(
    name='snn_pipeline',
    version='0.2.5',
    packages=find_packages(),
    package_data={'snn_pipeline': ['include data/*.pep', 'include data/*.bed', 'include data/*.hmm',
                                   'include data/species_list']},
    license='GPL-3.0',
    description='A pipeline for constructing and analyzing synteny network',
    long_description=open('README.md').read(),
    long_description_content_type='text/plain',
    author='Haochen li, Tao Zhao',
    author_email='lhc2018@nwafu.edu.cn, tao.zhao@nwafu.edu.cn',
    url='https://github.com/hcli007/SNN_pipeline',
    install_requires=['biopython>=1.79',
                      'igraph>=0.1.14',
                      'networkx>=2.6.3',
                      'numpy>=1.21.2',
                      'pandas>=1.4.1',
                      'scikit-learn>=1.3.2',
                      ],

    classifiers=['Development Status :: 4 - Beta',
                 'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
                 'Programming Language :: Python :: 3',
                 'Operating System :: POSIX :: Linux',
                 'Topic :: Scientific/Engineering :: Bio-Informatics',
                 ],
    python_requires='>=3.8',
    entry_points={
        'console_scripts': [
            'synetbuild = snn_pipeline.synetbuild:main',
            'synetfind = snn_pipeline.synetfind:main',
            'synetprefix = snn_pipeline.synetprefix:main',
            'synetcontext = snn_pipeline.synetcontext:main',
            'synetmapping = snn_pipeline.synetmapping:main'
        ],
    },
)
