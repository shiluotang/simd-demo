ACLOCAL_AMFLAGS = -I m4

AM_CPPFLAGS = -g -Wall
AM_CFLAGS	=
AM_CXXFLAGS	=
AM_LDFLAGS	=

noinst_PROGRAMS = simd-enabled simd-disabled

simd_enabled_CPPFLAGS = $(AM_CPPFLAGS) $(SIMD_FLAGS) -funroll-loops
simd_enabled_SOURCES = ./src/simd-enabled.cpp

simd_disabled_CPPFLAGS = $(AM_CPPFLAGS) -funroll-loops
simd_disabled_SOURCES = ./src/simd-disabled.cpp
