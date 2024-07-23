from setuptools import setup
from distutils.core import Extension
from setuptools import find_packages
from distutils.sysconfig import get_python_inc

import numpy as np

pymultilayer_module=Extension('_pymultilayer',
        sources=['multilayer.cpp',
#            'multilayer.h',
#            'pymultilayer.h',
            'pymultilayer.cpp',
            'pymultilayer_wrap.cxx',
            ],
        include_dirs=[np.get_include(),get_python_inc()],
        )

setup(name="pymultilayer",
        version="0.1",
        author="xxx",
        description="""pymultilayer module""",
        ext_modules=[pymultilayer_module],
        py_modules=["pymultilayer"],
        install_requires=['numpy'],
        packages=find_packages(include=["numpy"],),
        )
