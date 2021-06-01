set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

include_directories(
  ../
  ./
  ./field
  ./main
  ./gf_complete
)
#add_subdirectory(field)
add_subdirectory(gf_complete)
link_directories(.)
if (WIN32)
  add_executable(NodeRecoveryCmd 
    node_recover_cmd.cpp
    Node.h
	ConstructionDefs.h
	FieldDefs.h
	FieldElement.h
  )
else()
  add_executable(NodeRecoveryCmd 
    node_recover_cmd.cpp
    Node.h
	ConstructionDefs.h
	FieldDefs.h
	FieldElement.h
)
endif()


target_link_libraries(NodeRecoveryCmd gf_complete)

if(MSVC)
  set_target_properties(NodeRecoveryCmd PROPERTIES LINK_FLAGS " /ignore:4099")
  set_target_properties(NodeRecoveryCmd PROPERTIES VS_DEBUGGER_WORKING_DIRECTORY "$(OutDir)")
endif()

set_target_properties(NodeRecoveryCmd PROPERTIES OUTPUT_NAME "NodeRecoveryCmd")
# To link libraries from current directory on Linux.
set_target_properties(NodeRecoveryCmd PROPERTIES INSTALL_RPATH ".")


if (CMAKE_BUILD_TYPE=DEBUG)
	set(OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/Debug")
else()
	set(OUT_DIR "${CMAKE_CURRENT_BINARY_DIR}/Release")
endif()

set(PACKAGE_DIR ./)