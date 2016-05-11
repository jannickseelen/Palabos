# List all header files in palabos/src
# subdirs of palabos/src
set(sub_dirs "dataProcessors" "offLattice" "algorithm" "finiteDifference"
"atomicBlock" "io" "basicDynamics" "latticeBoltzmann"	"boundaryCondition" "libraryInterfaces"
"coProcessors" "multiBlock" "parallelism" "complexDynamics" "multiGrid" "particles" "core" "multiPhysics" "modules")
# List all header files in the sub directories
foreach(sub ${sub_dirs})
	file(GLOB FOUND_H_FILES "${sub}/*.h")
	if(FOUND_H_FILES)
		set(H_FILES ${H_FILES} ${FOUND_H_FILES})
	endif()
	file(GLOB FOUND_HH_FILES "${sub}/*.hh")
	if(FOUND_HH_FILES)
		set(HH_FILES ${HH_FILES} ${FOUND_HH_FILES})
	endif()
	file(GLOB FOUND_CXX_FILES "${sub}/*.cpp")
	if(FOUND_CXX_FILES)
		set(CXX_FILES ${CXX_FILES} ${FOUND_CXX_FILES})
	endif()
	set(DIRECTORY_LIST ${DIRECTORY_LIST} "${sub}/")
endforeach()

add_library(palabos SHARED ${CXX_FILES} ${H_FILES} ${HH_FILES})

set_target_properties(palabos PROPERTIES LINKER_LANGUAGE CXX)

target_include_directories(palabos PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
$<INSTALL_INTERFACE:"${CMAKE_INSTALL_PREFIX}/lib/">  # <prefix>/include/mylib
)

target_link_libraries(palabos tinyxml dlib)
