'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os

os.environ["CC"] = "/share/apps/local/bin/g++"
os.environ["CXX"] = "/share/apps/local/bin/g++"

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
smc_parallel = Extension(name="smc_parallel",sources=['./smc_parallel.cpp','dcdio.c', 'smc_core.cpp'],
                   include_dirs =
                         [numpy_include,'./','/share/apps/local/git_working_copies/test_open_acc/sphere_mc/extensions'],
                   library_dirs =
                         ['/share/apps/local/git_working_copies/test_open_acc/sphere_mc/extensions','/share/apps/local/pgi/linux86-64/16.10/lib/','/share/apps/local/pgi/16.10/share_objects/lib64/','/state/partition1/apps/local/pgi/linux86-64/16.10/lib/','/usr/lib/gcc/x86_64-redhat-linux/4.4.7/'],
                   libraries = ["oacc_smc","accapi", "accg", "accn", "accg2", "dl", "cudadevice", "pgmp", "numa", "pthread", "nspgc", "pgc", "m", "gcc", "c", "gcc"],
                   #libraries = ["oacc_smc"],
                   extra_compile_args = ["-std=c++11"]
                   )

# NumyTypemapTests setup
setup(  name        = "PR_EXTENSION",
        description = "Module is middle code between python and pgc++ openacc",
        author      = "Joseph E. Curtis",
        version     = "0.1",
        ext_modules = [smc_parallel]
        )

