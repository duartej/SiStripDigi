#############################################################################
#
# CMAKE build setup for SiStripDigi
#
# For building Marlin with cmake type:
# (1) $ mkdir build
# (2) $ cd build
# (3) $ cmake -C ../BuildSetup.cmake ..
# (4) $ make install
#
#
#############################################################################


#############################################################################
# Setup path variables
#############################################################################

# ILCSOFT base directory
SET( ILC_HOME "$ENV{ILCSOFT_HOME}"
    CACHE PATH "Path to ILC Software" FORCE )
    
# Path to CMake Modules
SET( CMAKE_MODULE_PATH "/usr/local/ilcsoft/v01-08/ilcutil/trunk/cmakemodules"
    CACHE PATH "Path to CMake Modules" FORCE )    

# Path to CLHEP, respectively to 'clhep-config'
SET( CLHEP_HOME "$ENV{CLHEP_BASE_DIR}"
     CACHE PATH "Path to CLHEP" FORCE )
     
# Path to GEAR
SET( GEAR_HOME "$ENV{GEAR}"
    CACHE PATH "Path to GEAR" FORCE )

# Path to GSL, respectively to 'gsl-config'                                  
#SET( GSL_HOME "/usr/bin"                                                    
#     CACHE PATH "Path to GSL" FORCE ) 
     
# Path to LCIO
SET( LCIO_HOME "$ENV{LCIO}"
    CACHE PATH "Path to LCIO" FORCE )

# Path to Marlin
SET( Marlin_HOME "$ENV{MARLIN}"
     CACHE PATH "Path to Marlin" FORCE )

# Path to ROOT, respectively to 'root-config'
SET( ROOT_HOME "/usr/local"
    CACHE PATH "Path to ROOT-CONFIG" FORCE )


#############################################################################
# Project Options
#############################################################################

# Project depends on ...
SET( PROJECT_DEPENDS "CLHEP GEAR LCIO Marlin ROOT" 
     CACHE STRING "SiStripDigi dependence" FORCE )
    
# Set CMAKE build type (None Debug Release RelWithDebInfo MinSizeRel), default: RelWithDebInfo
SET( CMAKE_BUILD_TYPE "RelWithDebInfo" CACHE STRING "Choose the type of build" FORCE )

# Set basic options
SET( INSTALL_DOC OFF CACHE BOOL "Set to ON not to skip build/install Documentation" FORCE )

# Advanced options
#SET( BUILD_SHARED_LIBS OFF CACHE BOOL "Set to OFF to build static libraries" FORCE )

