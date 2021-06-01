set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")


include_directories(
  ../
  ./
)
# CMake instructions to make the static lib
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
			 gf_wgen.c)