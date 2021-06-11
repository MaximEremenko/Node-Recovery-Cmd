if (MSVC)
add_compile_options(/Ox /DINTEL_SSE4 /DINTEL_PCLMUL /QxHost)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /Ox")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /O0")
else()
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
add_compile_options(-msse4 -DINTEL_SSE4 -maes -mpclmul -DINTEL_PCLMUL)
endif()


include_directories(
  ../
  ./
)

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