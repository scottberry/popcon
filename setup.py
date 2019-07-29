from distutils.core import setup, Extension
import pybind11

cpp_args = ['-std=c++11']

ext_modules = [
    Extension(
        'euclidean_distance',
        ['euclidean_distance.cpp'],
        include_dirs=[pybind11.get_include(),
                      pybind11.get_include(user=True)],
        libraries=['gmp'],
        language='c++',
        extra_compile_args=cpp_args
    ),
]

setup(
    name='popcon',
    version='0.0.1',
    author='Scott Berry',
    author_email='scottdberry@gmail.com',
    description='Implementation of euclidean distance in C',
    install_requires=[
        'pybind11>=2.2.1',
	'numpy',
	'pandas',
	'scipy'
    ],
    setup_requires=['pybind11>=2.2.1'],
    ext_modules=ext_modules,
)
