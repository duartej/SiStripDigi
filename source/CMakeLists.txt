#############################################################################
#
# CMAKE file for building SiStripDigi
#
#
#############################################################################
INCLUDE_DIRECTORIES( BEFORE "./include" )

# require proper c++
ADD_DEFINITIONS( "-Wall -ansi -pedantic" )
ADD_DEFINITIONS( "-Wno-long-long" )


# definitions
IF( CMAKE_BUILD_TYPE STREQUAL "Debug" OR
		CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo" )
	ADD_DEFINITIONS( "-DSISTRIPDIGIDEBUG" )
ENDIF()

# get list of all source files
AUX_SOURCE_DIRECTORY( ./src library_sources )

# include directories
INCLUDE_DIRECTORIES( ./include )
INSTALL_DIRECTORY( ./include DESTINATION include FILES_MATCHING PATTERN "*.h" )


ADD_SHARED_LIBRARY( SiStripDigi ${library_sources} )
INSTALL_SHARED_LIBRARY( SiStripDigi DESTINATION lib )
TARGET_LINK_LIBRARIES( SiStripDigi ${CMAKE_DL_LIBS} )

#DISPLAY_STD_VARIABLES()

