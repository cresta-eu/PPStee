##############################################################################
# cmake behaviour
##############################################################################

cmake_minimum_required(VERSION 2.6)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")


##############################################################################
# basic project settings
##############################################################################
project(ppstee)
set (ppstee_VERSION_MAJOR 0)
set (ppstee_VERSION_MINOR 3)
set (ppstee_VERSION_PATCH "0d")

set (PPSTEE_LINK_LIBRARIES)


##############################################################################
# project options
##############################################################################

option(PPSTEE_ENABLE_PARMETIS "Enable ParMETIS within PPStee" ON)
option(PPSTEE_ENABLE_PTSCOTCH "Enable PTScotch within PPStee" ON)
option(PPSTEE_ENABLE_ZOLTAN "Enable Zoltan within PPStee" ON)


##############################################################################
# setup MPI
##############################################################################
#set (CMAKE_C_COMPILER mpicc)
#set (CMAKE_CXX_COMPILER mpicxx)

find_package(MPI REQUIRED)
set (CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_CXX_COMPILE_FLAGS})
set (CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_CXX_LINK_FLAGS})
include_directories (${MPI_CXX_INCLUDE_PATH})
set (PPSTEE_LINK_LIBRARIES ${MPI_CXX_LIBRARIES} ${PPSTEE_LINK_LIBRARIES})


##############################################################################
# partitioners
##############################################################################

if(PPSTEE_ENABLE_PARMETIS)
  find_package(Parmetis REQUIRED)
  add_definitions(-DPPSTEE_WITH_PARMETIS)
  include_directories(${PARMETIS_INCLUDE_DIRS})
  set (PPSTEE_LINK_LIBRARIES ${PARMETIS_LIBRARIES} ${PPSTEE_LINK_LIBRARIES})
endif(PPSTEE_ENABLE_PARMETIS)

if(PPSTEE_ENABLE_PTSCOTCH)
  find_package(PTScotch REQUIRED)
  add_definitions(-DPPSTEE_WITH_PTSCOTCH)
  include_directories(${PTSCOTCH_INCLUDE_DIRS})
  set (PPSTEE_LINK_LIBRARIES ${PTSCOTCH_LIBRARIES} ${PPSTEE_LINK_LIBRARIES})
endif(PPSTEE_ENABLE_PTSCOTCH)

if(PPSTEE_ENABLE_ZOLTAN)
  find_package(Zoltan REQUIRED)
  add_definitions(-DPPSTEE_WITH_ZOLTAN)
  include_directories(${ZOLTAN_INCLUDE_DIRS})
  set (PPSTEE_LINK_LIBRARIES ${ZOLTAN_LIBRARIES} ${PPSTEE_LINK_LIBRARIES})
endif(PPSTEE_ENABLE_ZOLTAN)


##############################################################################
# uninstall target
##############################################################################
CONFIGURE_FILE(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake"
  IMMEDIATE @ONLY)
ADD_CUSTOM_TARGET(uninstall
  "${CMAKE_COMMAND}" -P "${CMAKE_CURRENT_BINARY_DIR}/cmake/cmake_uninstall.cmake") 


##############################################################################
# doc target
##############################################################################

find_package (Doxygen)
if (DOXYGEN_FOUND)
configure_file (
	"${ppstee_SOURCE_DIR}/doxygen.conf.in"
	"${ppstee_BINARY_DIR}/doxygen.conf"
	@ONLY
)
add_custom_target (doc
#	ALL					# adds doc to target ALL; uncomment if not wanted
	COMMAND ${DOXYGEN_EXECUTABLE} "${ppstee_BINARY_DIR}/doxygen.conf"
	COMMENT "Building doxygen documentation...")
endif(DOXYGEN_FOUND)


##############################################################################
# check-style target
##############################################################################

add_custom_target (check-style
	COMMAND ./audit.py xml
  COMMENT "Checking style...")


##############################################################################
# debug flag
##############################################################################

if (CMAKE_BUILD_TYPE MATCHES DEBUG)
  add_definitions(-DPPSTEE_DEBUG)
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)

option(GDBDEBUG "Build example with infinite loop at start to enable attaching gdb." "0")
add_definitions(-DGDBDEBUG=${GDBDEBUG})


##############################################################################
# main project settings
##############################################################################

# set some dirs
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${ppstee_SOURCE_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${ppstee_SOURCE_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${ppstee_SOURCE_DIR}/bin")

# external libs
#set (TRILIBS isorropia zoltan aztecoo triutils teuchos epetra ptscotch ptscotcherr ptscotcherrexit parmetis metis blas lapack)
#set (PETSCLIBS petsc zoltan ptscotch ptscotcherr ptscotcherrexit parmetis metis ml blas lapack X11)

# include
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/include")

# libraries (including their headers)
set (ppsteeLibSourceNames "ppstee" "ppstee_objects" "ppstee_error")

foreach (libname ${ppsteeLibSourceNames})
  set (ProjLibs ${ProjLibs} "src/${libname}.cpp")
  install (FILES "include/${libname}.hpp" DESTINATION include)
endforeach(libname)

add_library (ppstee ${ProjLibs})
install (TARGETS ppstee DESTINATION lib)

# examples
add_executable(ppstee_example "examples/ppstee_example.cpp")
target_link_libraries (ppstee_example ppstee ${PPSTEE_LINK_LIBRARIES})
install (TARGETS ppstee_example DESTINATION bin)

#add_executable (DlrMatrixTrilinosCrs_example "examples/dlrmatrixtrilinoscrs_example.cpp")
#target_link_libraries (DlrMatrixTrilinosCrs_example DlrMatrix ${TRILIBS} ${MPI_LIBRARIES})
#install (TARGETS DlrMatrixTrilinosCrs_example DESTINATION bin)


##############################################################################
# tests
##############################################################################

# set sources
set(TEST_SOURCES
  tools/googletest-1.6/gmock-gtest-all.cc
  test/mainTest.cpp
  test/MainSystemTest.cpp
  )
#set(TEST_SOURCES ${TEST_SOURCES} ${ProjLibs})

# include
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}/tools/googletest-1.6")
include_directories ("${CMAKE_CURRENT_SOURCE_DIR}")

# activate testing
enable_testing()

# define tests
if (CMAKE_BUILD_TYPE MATCHES DEBUG)
  set(TEST_MPI_NUMPROCS 3 113)
  #set(TEST_MPI_NUMPROCS 2 3 4 5 6)
else()
  set(TEST_MPI_NUMPROCS 2 3 4 5 6 12)
endif(CMAKE_BUILD_TYPE MATCHES DEBUG)
foreach(CURRENT_NUMPROC ${TEST_MPI_NUMPROCS})
  add_executable(ppstee_test_${CURRENT_NUMPROC} ${TEST_SOURCES})
  set_target_properties(ppstee_test_${CURRENT_NUMPROC} PROPERTIES COMPILE_DEFINITIONS USING_NUMPROCS_${CURRENT_NUMPROC})
  target_link_libraries(ppstee_test_${CURRENT_NUMPROC} ppstee ${PPSTEE_LINK_LIBRARIES})
  add_test(NAME MainSystemTest-${CURRENT_NUMPROC}
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${CURRENT_NUMPROC} $<TARGET_FILE:ppstee_test_${CURRENT_NUMPROC}>)
endforeach(CURRENT_NUMPROC)

# Adding custom test target "check" because it is the only way to show the full output of test errors on console
add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND} --force-new-ctest-process --output-on-failure)

##############################################################################
# CPack TGZ targets
##############################################################################

# configure package generator
set(CPACK_GENERATOR "TGZ")

# configure source package generator
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "ppstee-${ppstee_VERSION_MAJOR}.${ppstee_VERSION_MINOR}.${ppstee_VERSION_PATCH}")

#################
## set files/dirs/patterns to be included in source package
#set(SOURCE_PACKAGE_INCLUDE_FILES
#  "logo.png"
#  "PPStee_flow.png"
#  "CMakeLists.txt"
#  "/include/"
#  "/cmake/"
#  "/examples/"
#  "/tools/"
#  "/test/res/"
#  ${ProjLibs}
#  ${TEST_SOURCES}
#)
#list(REMOVE_DUPLICATES SOURCE_PACKAGE_INCLUDE_FILES)
#
## set files/dirs/patterns to be ignores in source package
#set(SOURCE_PACKAGE_IGNORE_FILES
#  "\\\\.exe"
#  "\\\\.tar\\\\.gz"
#)
#
## assemble CPACK_SOURCE_IGNORE_FILES
#foreach(SRCFILE ${SOURCE_PACKAGE_IGNORE_FILES})
#  list(APPEND SOURCE_PACKAGE_FILES_LIST ".*${SRCFILE}")
#endforeach(SRCFILE)
#foreach(SRCFILE ${SOURCE_PACKAGE_INCLUDE_FILES})
#  list(APPEND SOURCE_PACKAGE_FILES_LIST "^.*${SRCFILE}")
#endforeach(SRCFILE)
#list(APPEND SOURCE_PACKAGE_FILES_LIST ".*")
#set(CPACK_SOURCE_IGNORE_FILES "${SOURCE_PACKAGE_FILES_LIST}")
#################

set(CPACK_SOURCE_IGNORE_FILES "/CMakeFiles/;CMakeCache.txt;compile_commands.json$;cmake_install.cmake;.ycm_extra_conf.py.?;vera.*xml;cppcheck.xml;/cmake/cmake_uninstall.cmake$;install_manifest.txt;doxygen.conf$;Makefile;CPackConfig.cmake;/_CPack_Packages/;CPackSourceConfig.cmake;.*.exe$;.*.swp$;.*/.nfs.*;/Testing/;CTestTestfile.cmake;/*.*tar.gz;/inst/;/doc/;/bin/;/lib/;/.svn/;breakmain.gdb;debug-manually.sh;debug.sh;.cproject;.project;runlog;tags")

include(CPack)
