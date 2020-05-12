# From: https://github.com/hsidky/SAPHRON
#
# The MIT License (MIT)
#
# Copyright (c) 2015 Hythem Sidky
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Find QUADMATH.

SET(QUADMATH_NAMES ${QUADMATH_NAMES} libquadmath.a quadmath)
FIND_LIBRARY(QUADMATH_LIBRARY
  NAMES ${QUADMATH_NAMES}
  PATHS /usr/lib64/atlas /usr/lib/atlas 
        /usr/lib64 /usr/lib /usr/local/lib64 
        /usr/local/lib /usr/x86_64-linux-gnu/*
        /usr/lib/gcc/x86_64-linux-gnu/*
  )

IF (QUADMATH_LIBRARY)
  SET(QUADMATH_LIBRARIES ${QUADMATH_LIBRARY})
  SET(QUADMATH_FOUND "YES")
ELSE (QUADMATH_LIBRARY)
  SET(QUADMATH_FOUND "NO")
ENDIF (QUADMATH_LIBRARY)

IF (QUADMATH_FOUND)
   IF (NOT QUADMATH_FIND_QUIETLY)
      MESSAGE(STATUS "Found QUADMATH: ${QUADMATH_LIBRARIES}")
   ENDIF (NOT QUADMATH_FIND_QUIETLY)
ELSE (QUADMATH_FOUND)
   IF (QUADMATH_FIND_REQUIRED)
      MESSAGE(FATAL_ERROR "Could not find QuadMath")
   ENDIF (QUADMATH_FIND_REQUIRED)
ENDIF (QUADMATH_FOUND)