set(test_headers
	${CMAKE_CURRENT_SOURCE_DIR}/test_basic.h
	PARENT_SCOPE)

# Includ this directory for test headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
