from setuptools import setup

kwds = {'name': 'pyNURBS',
        'version': '1.0',
        'packages': ['pynurbs',
                     'pynurbs.geometry', 'pynurbs.geometry.methods',
                     'pynurbs.graphics', 'pynurbs.graphics.methods'],
        'install_requires': ['numpy', 'scipy', 'mayavi'],
        'author': 'trelau',
        'description': 'Simple Python-based NURBS library.',
        'license': 'BSD'}

setup(**kwds)
