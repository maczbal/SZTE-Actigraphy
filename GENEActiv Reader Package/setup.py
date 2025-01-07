from setuptools import setup, Extension
import pybind11

cpp_args = ['-std=c++11', '-stdlib=libc++', '-mmacosx-version-min=10.7']

sfc_module = Extension(
    'pybindRetTest',
    sources=['module.cpp'],
    include_dirs=[pybind11.get_include()],
    language='c++',
    extra_compile_args=cpp_args,
    )

setup(
    name='pybindRetTest',
    version='1.0',
    description='Python package with superfastcode C++ extension (PyBind11)',
    ext_modules=[sfc_module],
)