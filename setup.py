from setuptools import setup, find_packages
from scCirRL.__init__ import __program__, __module__, __version__, __author__, __email__, __description__, __url__
from scCirRL.__init__ import __scripts__


def read(path):
    with open(path, "r") as f:
        return f.read()
def read_list(path):
    return read(path).splitlines()

entrypoints = ["{}={}.main:main".format(__module__, __module__)]
entrypoints.extend(["{}={}.{}:main".format(script, __module__, script) for script in __scripts__])

long_description = read("README.md")
installrequires = read_list("pypi_requirements.txt")
setuprequires = ["setuptools", "wheel", "twine"]

setup(
    name=__program__,
    version=__version__,
    license="MIT",
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,
    description=__description__,
    long_description=long_description,
    long_description_content_type="text/markdown",
    zip_safe=False,
    packages=find_packages(),
    # packages=find_packages(include=[__module__, "{}.*".format(__module__)]),
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=installrequires,
    setup_requires=setuprequires,
    url=__url__,
    keywords=[__program__],
    entry_points={
        "console_scripts":entrypoints, 
    },
    classifiers=[
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3.8",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Utilities",
    ],
)
