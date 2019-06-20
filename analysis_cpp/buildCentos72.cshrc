#!/bin/csh

alias setprompt 'set prompt="${cwd}% "' 
setprompt # to set the initial prompt 
alias cd 'chdir \!* && setprompt' 


########################################################################
# PREFIX and non-system programs/libraries
########################################################################

### prefix area
setenv PREFIX /u/group/clas/builds/centos7/trunk

### non-system builds of programs and libraries
setenv GCC /apps/gcc/6.4.0
setenv ROOT /apps/root/6.12.06
setenv CERN /apps/cernlib/x86_64_rhel7/2005
setenv PYTHON /apps/python/2.7.12
setenv SCONS /apps/scons

setenv PLUTO_BUILD /work/clas/clasg12/mkunkel/PLUTO/pluto_v5.42
setenv PLUTO_SRC $PLUTO_BUILD/src
setenv PLUTO_LIB $PLUTO_BUILD/lib
setenv PLUTO_INC $PLUTO_BUILD/src
setenv PLUTO_PLUGINS $PLUTO_BUILD/plugins
 
setenv ROOTSYS $ROOT/root
########################################################################
# PATH
########################################################################

setenv PATH .:${PREFIX}/bin
setenv PATH ${PATH}:${PREFIX}/scripts

setenv PATH ${PATH}:${GCC}/bin
setenv PATH ${PATH}:${ROOT}/root/bin
setenv PATH ${PATH}:${PYTHON}/bin
setenv PATH ${PATH}:${SCONS}/bin
setenv PATH ${PATH}:${PREFIX}/build/bin
### standard system paths
setenv PATH ${PATH}:/site/bin:/apps/bin
setenv PATH ${PATH}:/usr/bin:/bin:/usr/sbin:/sbin

setenv PATH ${PATH}:./bin:./build/bin

########################################################################
# LD_LIBRARY_PATH
########################################################################

### run-time library loading path
setenv LD_LIBRARY_PATH .:${PREFIX}/lib

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${GCC}/lib64
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${ROOT}/root/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PYTHON}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${SCONS}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PREFIX}/build/lib

########################################################################
#  JAVA Environment
#  ########################################################################
module load java_1.8
#  use groovy
setenv PATH /work/clas12/devita/groovy/2.4.7/bin/:$PATH
#
#  ########################################################################
#  # PATHS and Executables
#  ########################################################################
setenv CCDB /group/clas12/src/ccdb_0.08
setenv PATH ${PATH}:${CCDB}/bin
setenv CCDB_USER $USER
setenv CCDB_CONNECTION mysql://clas12reader@clasdb.jlab.org/clas12
setenv CCDBINC ${CCDB}/include
setenv CCDBLIB ${CCDB}/lib

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${CCDB}/lib
setenv PYTHONPATH ${CCDB}/python:${CCDB}/python/ccdb/ccdb_pyllapi/


setenv CLARA_HOME /work/clas12/segarrae/clara
setenv CLARA_USER_DATA /work/clas12/segarrae/clara_user_data
setenv COATJAVA $CLARA_HOME/plugins/clas12
setenv CLAS12DIR $CLARA_HOME/plugins/clas12
setenv PATH ${COATJAVA}/bin:${PATH}
#
setenv MAVEN /apps/maven/PRO
setenv PATH ${MAVEN}/bin:${PATH}

########################################################################
# PYTHONPATH
########################################################################

### python modules search path

#setenv PYTHONPATH ${PREFIX}/lib/python
#

########################################################################
# sources for build directories
########################################################################

setenv MYSQLINC /usr/include/mysql
setenv MYSQLLIB /usr/lib64/mysql

setenv TCLLIB /usr/lib64

setenv BOOSTINC /usr/include/boost
setenv BOOSTLIB /usr/lib64

setenv CERNLIB ${CERN}/lib




########################################################################
# misc
########################################################################

alias		clara		${CLARA_HOME}/bin/clara-shell

rehash

