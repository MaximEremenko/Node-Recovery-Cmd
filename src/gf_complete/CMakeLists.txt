#set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -msse4 -DINTEL_SSE4 -maes -mpclmul -DINTEL_PCLMUL")

if (MSVC)
    # warning level 4 and all warnings as errors
    #add_compile_options(/W4 /WX)
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox /DINTEL_SSE4 /DINTEL_PCLMUL /QxHost")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /O0 /DINTEL_SSE4 /DINTEL_PCLMUL /QxHost")
else()
    # lots of warnings and all warnings as errors
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -maes -mpclmul -msse4 -DINTEL_SSE4 -DINTEL_PCLMUL")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} -O0 -maes -mpclmul -msse4 -DINTEL_SSE4 -DINTEL_PCLMUL")
    #add_compile_options(-msse4 -DINTEL_SSE4 -maes -mpclmul -DINTEL_PCLMUL)
endif()


include_directories(
  ../
  ./
)
#if (WIN32)
#include_directories(
#  ./sys
#)
#endif()

# CMake instructions to make the static lib
if (MSVC)
add_library( gf_complete STATIC
             gf.c 
			 gf_complete.h
			 gf_general.h
			 gf_int.h
			 gf_method.c
			 gf_method.h
			 gf_rand.c
			 gf_rand.h
			 gf_w4.c
			 gf_w8.c
			 gf_w16.c
			 gf_w32.c
			 gf_w64.c
			 gf_w128.c
			 gf_wgen.c
			 ./sys/time.cpp
			 ./sys/time.h
			 ./sys/times.h
			 )
else()
add_library( gf_complete STATIC
             gf.c 
			 gf_complete.h
			 gf_general.h
			 gf_int.h
			 gf_method.c
			 gf_method.h
			 gf_rand.c
			 gf_rand.h
			 gf_w4.c
			 gf_w8.c
			 gf_w16.c
			 gf_w32.c
			 gf_w64.c
			 gf_w128.c
			 gf_wgen.c
			 )
endif()