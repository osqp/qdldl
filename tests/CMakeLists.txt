# Add subdirectories of all the tests
add_subdirectory(basic)

set(test_headers ${headers} PARENT_SCOPE)

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
