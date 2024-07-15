# Makefile
# copy general purpose .h and .c files here for distribution
# anything in this directory is just a copy from elsewhere (except Makefile)

all:
	@echo 'choose "get"'

#
# jf: copy files from master source directories to these local copies
#
home = $(HOME)
jf_def = $(home)/l/src/defs/
jf_thr = $(home)/l/src/util/thread/
jf_umx = $(home)/l/src/util/mex/

get:
	cp -pf \
		$(jf_def)/defs-env.h \
	\
		$(jf_thr)/jf,thread1.c \
		$(jf_thr)/jf,thread1.h \
	\
		$(jf_umx)/jf,mex,def.h \
		$(jf_umx)/mexarg.c \
		$(jf_umx)/def,mexarg.h .
