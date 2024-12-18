#######################################
# debian.mk
# Default options for debian computer
#######################################
CC=gcc
LIBSLOCAL=-L/usr/lib -llapack -lblas -lm
INCLUDEBLASLOCAL=-I/usr/include
OPTCLOCAL=-fPIC -march=native
