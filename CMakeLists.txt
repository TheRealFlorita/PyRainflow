cmake_minimum_required(VERSION 3.16 FATAL_ERROR)

set(PROJECT_NAME PyRainflow)
project(${PROJECT_NAME} LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(${CMAKE_INSTALL_PREFIX} "${CMAKE_CURRENT_SOURCE_DIR}/install/")

# Header and source files
set(Headers "Rainflow.h" "RamerDouglasPeucker.h")
set(Sources "Rainflow.cpp" "RamerDouglasPeucker.cpp" "Main.cpp")

# Target
add_library(${PROJECT_NAME} SHARED ${Headers} ${Sources})

if(MSVC)
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "/O2 /arch:AVX /EHsc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /O2 /arch:AVX /EHsc")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /O2 /arch:AVX /EHsc")
else()
    set_target_properties(${PROJECT_NAME} PROPERTIES COMPILE_FLAGS "-O3 -mavx")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -mavx")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -mavx")
endif()

# Find Python3 package
find_package(Python3 REQUIRED COMPONENTS Development)
include_directories(${Python3_INCLUDE_DIRS})

# Link Python3 library
target_link_libraries(${PROJECT_NAME} PUBLIC ${Python3_LIBRARIES})
target_link_directories(${PROJECT_NAME} PUBLIC ${Python3_LIBRARY_DIRS})

# Install target
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})