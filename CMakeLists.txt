cmake_minimum_required(VERSION 3.24)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
   message(FATAL_ERROR "In-source builds are not allowed!")
endif()

set(TARGET_NAME BPOTF)

# Project name and options
project(
   ${TARGET_NAME}
   LANGUAGES CXX
)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")  # or 11, 14, 17, 20
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # ensure standard is supported
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
   message(INFO " No optimization!")
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -Wall -Wstack-protector -flto=auto -DDEBUG_OBPOTF")
else()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wstack-protector -flto=auto")
endif()

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_RUNTIME_OUTPUT_DIRECTORY})
set(BUILD_OBJECTS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/.obj)

include(FetchContent)
FetchContent_Declare(
   pybind11
   GIT_REPOSITORY https://github.com/pybind/pybind11
   GIT_TAG        master
   OVERRIDE_FIND_PACKAGE
)
find_package(pybind11 REQUIRED)

include_directories(
   ${CMAKE_CURRENT_LIST_DIR}/src
)

file(
   GLOB_RECURSE
   BPOTF_SRCS
   CMAKE_CONFIGURE_DEPENDS
   "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp"
)

list(REMOVE_ITEM BPOTF_SRCS "${CMAKE_CURRENT_LIST_DIR}/src/CSC/main_test.cpp")
list(REMOVE_ITEM BPOTF_SRCS "${CMAKE_CURRENT_LIST_DIR}/src/CSR/main_test.cpp")

pybind11_add_module(${TARGET_NAME} ${BPOTF_SRCS})

add_custom_command(
   TARGET ${TARGET_NAME} POST_BUILD
   COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${TARGET_NAME}> ${CMAKE_CURRENT_LIST_DIR}/src
)