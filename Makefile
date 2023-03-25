FC=gfortran
FFLAGS=-O3 -Wall -Wextra -std=f2008
CONDAPATH="False"
OUTFILE=fortplanet

PROJDIR := $(realpath $(CURDIR)/..)
SRCDIR := ./lib/fortran
TARGETDIR := ./lib/fortran

_SRC=my_types.f95 run_params.f95 constants.f95 LinAlg.f95 class_table_new.f95 phase.f95 eosfort.f95 utils.f95 integrator.f95 eosmat.f95 fortshell.f95 fortlayer.f95 fortplanet.f95
_OBJ=${_SRC:.f95=.o}

SRC=$(patsubst %, ${SRCDIR}/%,${_SRC})
OBJ=$(patsubst %, ${TARGETDIR}/%,${_OBJ})

printvars:
	@echo "PROJDIR:" ${PROJDIR}
	@echo "SRCDIR:" ${SRCDIR}
	@echo "SRC:" ${SRC}
	@echo "OBJ:" ${OBJ}
    
%.o: %.f95
	${FC} ${FFLAGS} -o $@ -c $<

#Create static library from source files
static:
	gfortran -g -fbacktrace -fPIC -c ${SRC}
	@mv ./*.o ./*.mod  ${TARGETDIR}

#Create python wrapper from static library using f2py
wrapper:
	python3 -m numpy.f2py -c -I${TARGETDIR} ${SRCDIR}/eosfort_wrapper.f95 ${OBJ} -m ${OUTFILE}
	f2py -c -I${TARGETDIR} ${SRCDIR}/functionsPy.f95 ${OBJ} -m fortfunctions
	@mkdir -p ${TARGETDIR}
	@mv ./*.so ${TARGETDIR}
	
clean:
	@rm -f ${SRCDIR}/*.mod ${SRCDIR}/*.o

install: static wrapper clean
	@if(${CONDAPATH}="True");then\
		conda develop ${CURDIR}/lib;\
		conda develop ${CURDIR}/src/python;\
	fi


