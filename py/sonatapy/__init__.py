# get the version
from ._version import __version__

# explicitly set the package variable to ensure relative import work
__package__ = "sonatapy"

# import important stuff
from .bok import Bok
