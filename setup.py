from setuptools import setup

setup(
        name='conformer',
        version='0.91',
        py_modules=['conformer'],
        install_requires=[
            'datetime', 'argparse', 'rdkit', 
            ],
        entry_points='''
        [console_scripts]
        conformer = conformer:conformer
        ''',
        )
