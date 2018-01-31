add_subdirectory(external/src/neb)
add_library(cppneb cppneb_wrapper.f90 neb.f90)
target_link_libraries(cppneb neb)
