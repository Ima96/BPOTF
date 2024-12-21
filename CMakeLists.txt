cmake_minimum_required(VERSION 3.24)

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
   message(FATAL_ERROR "In-source builds are not allowed!")
endif()

set(TARGET_NAME BPOTF)

if(NOT DEFINED CPP_BPOTF_VER)
   set(CPP_BPOTF_VER "0.0.dev0")
endif(NOT DEFINED CPP_BPOTF_VER)

string(REGEX MATCH "([0-9]*)\.([0-9]*)\.([0-9a-z]*)" _ ${CPP_BPOTF_VER})
set(version_major ${CMAKE_MATCH_1})

# string(REGEX MATCH "VERSION_MINOR ([0-9]*)" _ ${ver})
set(version_minor ${CMAKE_MATCH_2})

# string(REGEX MATCH "VERSION_PATCH ([0-9]*)" _ ${ver})
set(version_patch ${CMAKE_MATCH_3})

# Project name and options
project(
   ${TARGET_NAME}
   LANGUAGES CXX
   VERSION ${version_major}.${version_minor}
)

set(CMAKE_CXX_STANDARD 20 CACHE STRING "C++ version selection")  # or 11, 14, 17, 20
set(CMAKE_CXX_STANDARD_REQUIRED ON)  # ensure standard is supported

if (MSVC)
   if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      add_compile_options(/Od /RTC1 /Zi /MDd)
   else()
      # WX warning level: 3 equivalent to -Wall, 4 -Wall -Wextra
      # GS is for stack protection like -Wstack-protector
      add_compile_options(/O2 /GL /Ob2 /W3 /GS /MD)
      # LTCG enables link time Code generation for improvements
      add_link_options(/LTCG)
   endif()

   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
   set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

else()
   if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      message(INFO " No optimization!")
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -Wall -Wstack-protector -flto=auto -DDEBUG_OBPOTF")
   else()
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -Wstack-protector -flto=auto")
   endif()

endif()

option(DEBUG_TIMING  "Enable timing prints." OFF)

if(DEBUG_TIMING)
   add_definitions(-DTIMING_OBPOTF)
endif()

if(NOT SETUP_PY_BUILD)
   set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${TARGET_NAME})
endif()
set(BUILD_OBJECTS_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/.obj)

include(FetchContent)
FetchContent_Declare(
   pybind11
   GIT_REPOSITORY https://github.com/pybind/pybind11
   GIT_TAG        master
   GIT_SHALLOW    TRUE
)
FetchContent_MakeAvailable(pybind11)

include_directories(
   ${CMAKE_CURRENT_LIST_DIR}/src
)

# Gather sources
file(
   GLOB_RECURSE
   BPOTF_SRCS
   CMAKE_CONFIGURE_DEPENDS
   "${CMAKE_CURRENT_LIST_DIR}/src/*.cpp"
)

# Remove test files
list(REMOVE_ITEM BPOTF_SRCS "${CMAKE_CURRENT_LIST_DIR}/src/CSC/main_test.cpp")
list(REMOVE_ITEM BPOTF_SRCS "${CMAKE_CURRENT_LIST_DIR}/src/CSR/main_test.cpp")

# Try to detect automatic version change
set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_CURRENT_LIST_DIR}/version)

# Generate the module
pybind11_add_module(${TARGET_NAME} ${BPOTF_SRCS})
message(NOTICE "CMAKE_LIBRARY_OUTPUT_DIRECTORY = ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

# Add compile definitions
target_compile_definitions(${TARGET_NAME} PUBLIC BPOTF_VERSION="${CPP_BPOTF_VER}")

if (MSVC)
   # Ensure the target output directory is set for all configurations for MSVC because
   # it generates multiple configurations and it places the output on another folder.
   set_target_properties(${TARGET_NAME} PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
      RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
      RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}"
   )
endif()

if (NOT SETUP_PY_BUILD OR EDITABLE_INST)
   add_custom_command(
      TARGET ${TARGET_NAME} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${TARGET_NAME}> ${CMAKE_SOURCE_DIR}/src
      COMMENT "Copying module $<TARGET_FILE:${TARGET_NAME}> to ${CMAKE_SOURCE_DIR}/src"
   )
endif()