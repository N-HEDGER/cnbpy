from os.path import join as pjoin

# Format expected by setup.py and doc/source/conf.py: string of form "X.Y.Z"
_version_major = 0
_version_minor = 1
_version_micro = ''  # use '' for first of series, number for 1 and above
_version_extra = 'dev'
# _version_extra = ''  # Uncomment this for full releases

# Construct full version string from these.
_ver = [_version_major, _version_minor]
if _version_micro:
    _ver.append(_version_micro)
if _version_extra:
    _ver.append(_version_extra)

__version__ = '.'.join(map(str, _ver))

CLASSIFIERS = ["Development Status :: 3 - Alpha",
               "Environment :: Console",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: GPL-v3 License",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Topic :: Scientific/Engineering"]

# Description should be a one-liner:
description = "cnbpy: a package for interacting with ABIDE 2 data"

# Long description will go up on the pypi page
long_description = """
CNBPY
========
Placeholder
License
=======
Placeholder
"""

NAME = "CNBPY"
MAINTAINER = "Nicholas Hedger"
MAINTAINER_EMAIL = "nhedger1@gmail.com"
DESCRIPTION = description
LONG_DESCRIPTION = long_description
URL = "https://github.com/N-HEDGER/cnbpy"
DOWNLOAD_URL = ""
LICENSE = "GPL3"
AUTHOR = "Nicholas Hedger"
AUTHOR_EMAIL = "nhedger1@gmail.com"
PLATFORMS = "OS Independent"
MAJOR = _version_major
MINOR = _version_minor
MICRO = _version_micro
VERSION = __version__
PACKAGE_DATA = {'cnbpy': [pjoin('test', 'data', '*')]}
REQUIRES = ["numpy","scipy","pandas","datalad","nibabel"]
