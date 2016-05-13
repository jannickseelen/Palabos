#######################################
# CMAKE PALABOS MAIN
#######################################

# Force CMake version 3.1 or above
cmake_minimum_required (VERSION 2.8)

# Project name:
project(PALABOS)

# Palabos ROOT
set(PALABOS_ROOT ${CMAKE_CURRENT_SOURCE_DIR})

# Set Location of CMake scripts
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMake/")

#Add tinyxml from externalLibraries
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/tinyxml")
file(GLOB CXX_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/tinyxml/*.cpp")
file(GLOB C_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/tinyxml/*.c")
file(GLOB HH_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/tinyxml/*.hh")
file(GLOB H_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/tinyxml/*.h")
set(TINYXML_SOURCE ${CXX_FILES} ${C_FILES} ${H_FILES} ${HH_FILES})
add_library(tinyxml STATIC ${TINYXML_SOURCE})


#include X11
#find_package(X11 REQUIRED)
#include_directories(${X11_INCLUDE_DIRS})
#add_library(X11lib STATIC ${X11_LIBRARIES} ${X11_INCLUDE_DIRS})
#set_target_properties(X11lib PROPERTIES LINKER_LANGUAGE CXX)
#set(X11_INCLUDE_DIR "/usr/lib/x86_64-linux-gnu/")
#include_directories("${CMAKE_CURRENT_SOURCE_DIR}/externalLibraries/X11")
#file(GLOB CXX_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/X11/*.cpp")
#file(GLOB C_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/X11/*.c")
#file(GLOB HH_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/X11/*.hh")
#file(GLOB H_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/X11/*.h")
#set(X11_SOURCE ${CXX_FILES} ${C_FILES} ${H_FILES} ${HH_FILES})

# include dlib
#add_definitions("-DDLIB_ENABLE_ASSERTS")
#add_definitions("-DDLIB_ENABLE_STACK_TRACE")
#add_definitions("-DDLIB_ISO_CPP_ONLY")
#add_definitions("-DDLIB_NO_GUI_SUPPORT")
#set(dlib_threads "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/dlib/threads")
#include_directories(${dlib_threads})
#include("${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/dlib/cmake")
#file(GLOB CXX_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/dlib/all/source.cpp")
#file(GLOB HH_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/dlib/all/*.hh")
#file(GLOB H_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/externalLibraries/dlib/all/*.h")
#set(DLIB_SOURCE ${CXX_FILES} ${C_FILES} ${H_FILES} ${HH_FILES})
#add_library(dlib STATIC ${DLIB_SOURCE} ${X11_SOURCE})
#-DDLIB_ENABLE_ASSERTS=ON -DDLIB_ENABLE_STACK_TRACE=ON -DDLIB_ISO_CPP_ONLY=ON
#add_library(dlib SHARED ${DLIB_SOURCE})

#====================================================================================================================
#
#		OPTIONS
#
#====================================================================================================================


#SET(CLANG_LIB "/Library/Developer/CommandLineTools/usr/lib/")
#SET(CLANG_INCLUDE "/Library/Developer/CommandLineTools/usr/include/")

#SET(GNU_LIB "/usr/local/Cellar/gcc/5.3.0/lib/")
#SET(GNU_INCLUDE "/usr/local/Cellar/gcc/5.3.0/include/")

IF(CMAKE_COMPILER_IS_GNUCXX)
	SET(COMPILER_FLAGS "-std=gnu++14")
ELSE(CMAKE_COMPILER_IS_GNUCXX)
	SET(COMPILER_FLAGS "-std=c++11 -stdlib=libc++ -I${CLANG_INCLUDE} -L${CLANG_LIB}")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

SET(COMPILER_FLAGS "${COMPILER_FLAGS} -g -Wall")
SET(COMPILER_FLAGS "${COMPILER_FLAGS} -m64 -Wno-unused-variable -Wno-sign-compare")
SET(COMPILER_FLAGS "${COMPILER_FLAGS} -O0 -fdiagnostics-show-option")
SET(COMPILER_FLAGS "${COMPILER_FLAGS} -Wno-unused-value -Wno-deprecated")
# system dlib -L/usr/local/opt/openblas/lib  -I/usr/local/opt/openblas/include
# -O0 = optimization is 0 (debug mode) if -O4 (max optimized)
# -g = debug information
# -ferror-limit = 1 -> Max errors shown
# -arch x86_64 = force 64 bit architecture
# -DPLB_MAC_OS_X = Palabos Mac OSX flag
# -Wshorten-64-to-32 = Output warnings when 64bit code is truncated into 32 bit code
#  Force C++11 --> For CLANG use -std=c++11 -stdlib=libc++ -nostdinc++ FOR GNU use -std=gnu++11 -stdlib=libstdc++
# -Wimplicit-function-declaration = Warns if a function is used before it is declared (Linking error)
# -Wall = Turns on a lot of usefull warnings
# -ftemplate-backtrace-limit=0 See all backtraces
# -fno-common = Compile common globals like normal definitions
# -fobjc-gc = Enable Objective-C garbage collection
# -mthread-model <value> = The thread model to use, e.g. posix, single (posix by default)
# -pthread = Support POSIX threads in generated code

SET(DEFAULT ON CACHE INTERNAL "Default value for enabled by default options")
OPTION(ENABLE_MPI "Enable MPI" ${DEFAULT})

IF(APPLE)
	ADD_DEFINITIONS("-DPLB_MAC_OS_X")
	set(COMPILER_FLAGS "-DPLB_MAC_OS_X ${COMPILER_FLAGS}")
	message("THIS A MAC OSX SYSTEM... DEFINING PLB_MAC_OSX")
ENDIF(APPLE)

IF(ENABLE_MPI)
	INCLUDE(FindMPI)
	IF(MPI_C_FOUND AND MPI_CXX_FOUND)
		SET(CMAKE_C_COMPILER ${MPI_C_COMPILER})
		SET(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
		ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
		set(COMPILER_FLAGS "-DPLB_MPI_PARALLEL ${COMPILER_FLAGS}")
		message("MPI PARALLEL IS ENABLED... DEFINING PLB_MPI_PARALLEL")
		include_directories(${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
		add_library(openmpi STATIC ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
		set_target_properties(openmpi PROPERTIES LINKER_LANGUAGE CXX)
	ELSE(MPI_C_FOUND AND MPI_CXX_FOUND)
		#Check old FindMPI version
		IF(MPI_COMPILER)
			SET(CMAKE_C_COMPILER ${MPI_COMPILER})
			SET(CMAKE_CXX_COMPILER ${MPI_COMPILER})
			ADD_DEFINITIONS("-DPLB_MPI_PARALLEL")
			set(COMPILER_FLAGS "-DPLB_MPI_PARALLEL ${COMPILER_FLAGS}")
			message("MPI PARALLEL IS ENABLED... DEFINING PLB_MPI_PARALLEL")
			include_directories(${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
			add_library(openmpi STATIC ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
			set_target_properties(openmpi PROPERTIES LINKER_LANGUAGE CXX)
		ELSE(MPI_COMPILER)
			MESSAGE(FATAL_ERROR "MPI-COMPILER NOT found!")
		ENDIF(MPI_COMPILER)
	ENDIF(MPI_C_FOUND AND MPI_CXX_FOUND)
ENDIF(ENABLE_MPI)

OPTION(ENABLE_POSIX "Enable POSIX" ${DEFAULT})

IF(ENABLE_POSIX)
	ADD_DEFINITIONS("-DPLB_USE_POSIX ")
	ADD_DEFINITIONS("-pthread ")
	set(COMPILER_FLAGS "-DPLB_USE_POSIX -pthread ${COMPILER_FLAGS}")
	message("POSIX IS ENABLED... DEFINING PLB_USE_POSIX")
ENDIF(ENABLE_POSIX)

#OPTION(ENABLE_SMP_PARALLEL "Enable SMP_PARALLEL" ${DEFAULT})

IF(ENABLE_SMP_PARALLEL)
	ADD_DEFINITIONS("-DPLB_SMP_PARALLEL")
	ADD_DEFINITIONS("-mthread-model SMP ")
	set(COMPILER_FLAGS "-DPLB_SMP_PARALLEL -mthread-model SMP ${COMPILER_FLAGS}")
	message("SMP PARALLEL IS ENABLED... DEFINING PLB_SMP_PARALLEL")
ENDIF(ENABLE_SMP_PARALLEL)

#Set Debug flags
ADD_DEFINITIONS("-DPLB_DEBUG")
set(COMPILER_FLAGS "-DPLB_DEBUG  ${COMPILER_FLAGS}")

message("COMPILER_FLAGS: " ${COMPILER_FLAGS})
SET(CMAKE_BUILD_TYPE "Debug") #None, Debug, Release
SET(CMAKE_CXX_FLAGS ${COMPILER_FLAGS})
SET(CMAKE_C_FLAGS ${COMPILER_FLAGS})

#====================================================================================================================
#
#		FILES
#
#====================================================================================================================

# Select the executable cpp
set(PROJECT_PATH "/viscosityTest/src/")
set(PROJECT_EXE "test.cpp")

# Make a list for the source files
set(SOURCE "${PALABOS_ROOT}${PROJECT_PATH}${PROJECT_EXE}")
# Make a list for the include directories
set(DIRECTORY_LIST "${PALABOS_ROOT}${PROJECT_PATH}")

#Include the Tinyxml include dir
set(DIRECTORY_LIST ${DIRECTORY_LIST} ${TINYXML_INCLUDE_DIR})


add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/src/")



##################################################################################
#
#		LINK AND INSTALL
#
##################################################################################

add_executable(viscosityTest ${SOURCE})

target_include_directories(viscosityTest PUBLIC
$<BUILD_INTERFACE:${PROJECT_PATH}>
$<INSTALL_INTERFACE:"${CMAKE_INSTALL_PREFIX}/bin/">  # <prefix>/include/mylib
)


target_link_libraries(viscosityTest palabos tinyxml openmpi)
