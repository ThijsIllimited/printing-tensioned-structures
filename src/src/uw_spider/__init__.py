"""
********************************************************************************
uw_spider
********************************************************************************

.. currentmodule:: uw_spider


.. toctree::
    :maxdepth: 1


"""

from __future__ import print_function

import os


__author__ = ["Thijs Masmeijer, Tomas Mendez Echengucia, Ed Habtour"]
__copyright__ = "University of Washington"
__license__ = "MIT License"
__email__ = "thijsmasmeijer96@gmail.com, tmendeze@uw.edu, habtour@uw.edu "
__version__ = "0.1.0"


HERE = os.path.dirname(__file__)

HOME = os.path.abspath(os.path.join(HERE, "../../"))
DATA = os.path.abspath(os.path.join(HOME, "data"))
DOCS = os.path.abspath(os.path.join(HOME, "docs"))
TEMP = os.path.abspath(os.path.join(HOME, "temp"))


__all__ = ["HOME", "DATA", "DOCS", "TEMP"]
