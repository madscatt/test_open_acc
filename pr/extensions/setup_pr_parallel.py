'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os

#os.environ["CC"] = "/opt/pgi/osx86-64/16.10/bin/pgc++"
#os.environ["CXX"] = "/opt/pgi/osx86-64/16.10/bin/pgc++"
#os.environ["OPT"] = ""
#os.environ["CFLAGS"] = ""
#os.environ["CCSHARED"] = ""
#os.environ["LDSHARED"] = ""
#os.environ["SHLIB_SUFFIX"] = ""
#os.environ["AR"] = ""
#os.environ["ARFLAGS"] = ""

#os.environ["CFLAGS"] = "-c -static -Minfo=accel -ta=multicore -O3"

# System imports
from distutils.core import *
from distutils      import sysconfig

# Third-party modules - we depend on numpy for everything
import numpy

from numpy.distutils.core import Extension, setup

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

# simple extension module
pr_parallel = Extension(name="pr_parallel",sources=['./pr_parallel.c'],
                   include_dirs = [numpy_include,'./'],
                   #extra_link_args=['-framework', 'OpenGL', '-framework', 'GLUT'])
                   extra_objects = ["oacc_pr.a"],
                   )

# NumyTypemapTests setup
setup(  name        = "PR_EXTENSION",
        description = "Module is middle code between python and pgc++ openacc",
        author      = "Joseph E. Curtis",
        version     = "0.1",
        ext_modules = [pr_parallel]
        )

