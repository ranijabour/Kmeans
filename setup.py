
from distutils.core import setup, Extension
module = Extension("spkmeansmod",
                  sources=[
                    'spkmeansmodule.c',
                    'spkmeans.c'
                  ])
setup(
    name='spkmeansmod',
    author='Lama and Rani',
    version='1.0',
    description='kmeanscalc',
    ext_modules=[module])
