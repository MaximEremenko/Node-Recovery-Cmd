cmake_minimum_required (VERSION 3.9.0 FATAL_ERROR)

if (WIN32)
  set(ASM_LANG "ASM_MASM")
else()
  set(ASM_LANG "ASM")
endif()


project(NodeRecoveryCmd LANGUAGES CXX C ${ASM_LANG})

# Static linking
if (NOT WIN32)
  #SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
  SET(BUILD_SHARED_LIBS OFF)
  SET(BUILD_SHARED_LIBRARIES OFF)
  SET(CMAKE_EXE_LINKER_FLAGS "-static-libgcc -static-libstdc++")
endif()

# C++ flags and standarts
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CUDA_STANDARD 17)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)

# find_package(Threads REQUIRED)

set(CMAKE_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/bin/")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -stdlib=libc++ -I/usr/include/libcxxabi")
endif()

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

if(WIN32)
  add_definitions(-DNOMINMAX)
  add_definitions(-DUNICODE -D_UNICODE)
  add_definitions(-D_WIN32_WINNT=0x0A00)
endif()

string(TIMESTAMP BUILD_VERSION "%Y.%m.%d.%H" UTC)
execute_process(
  COMMAND git log -1 --format=%h
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  OUTPUT_VARIABLE GIT_COMMIT_HASH
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
string(CONCAT BUILD_VERSION ${BUILD_VERSION} "-" ${GIT_COMMIT_HASH})

# add all other projects
add_subdirectory(src)

set(PLD_PACKAGE_NAME "NodeRecoveryCmd-${BUILD_VERSION}")
#CPACK

set(CPACK_PACKAGE_VERSION ${BUILD_VERSION})
set(CPACK_PACKAGE_FILE_NAME ${PLD_PACKAGE_NAME})

if (WIN32)
	set(CPACK_GENERATOR "ZIP")
	set(CPACK_ZIP_COMPONENT_INSTALL ON)
else()
	set(CPACK_GENERATOR "TGZ")
	set(CPACK_TGZ_COMPONENT_INSTALL ON)
endif()

include (CPack)
