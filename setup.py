from setuptools import setup, find_packages

setup(
    name="sc_simulator",
    version="3.2",
    packages=find_packages(),
    py_modules=["__init__"],   # 让根目录直接变成模块
    # install_requires=[],
)