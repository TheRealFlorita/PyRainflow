cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

set(PROJECT_NAME PyRainflow)
project(${PROJECT_NAME} LANGUAGES C CXX)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(${CMAKE_INSTALL_PREFIX} "${CMAKE_CURRENT_SOURCE_DIR}/install/")

# Source groups
set(Headers "PyRainflow.h")
set(Sources "PyRainflow.cpp" "Source.def")

# Target
add_library(${PROJECT_NAME} SHARED ${Headers} ${Sources})
set_target_properties(${PROJECT_NAME} PROPERTIES VS_GLOBAL_KEYWORD "Win32Proj")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/Python3/include")

# Compile definitions
target_compile_definitions(${PROJECT_NAME} PRIVATE
	"$<$<CONFIG:Debug>:"
		"_DEBUG"
	">"
	"$<$<CONFIG:Release>:"
		"NDEBUG"
	">"
	"RAINFLOW_EXPORTS;"
	"_WINDOWS;"
	"_USRDLL;"
	"UNICODE;"
	"_UNICODE"
)

# Compile and link options
target_compile_options(${PROJECT_NAME} PRIVATE
	$<$<CONFIG:Debug>:
		/permissive-;
		/Od;
		/sdl;
		/W3
	>
	$<$<CONFIG:Release>:
		/std:c++17;
		/O2;
		/Ot
	>
        /arch:AVX;
	/EHsc;
	${DEFAULT_CXX_DEBUG_INFORMATION_FORMAT};
	${DEFAULT_CXX_EXCEPTION_HANDLING}
)

target_link_options(${PROJECT_NAME} PRIVATE
	$<$<CONFIG:Debug>:
		/INCREMENTAL
	>
	$<$<CONFIG:Release>:
		/OPT:REF;
		/OPT:ICF;
		/INCREMENTAL:NO
	>
	/DEBUG;
	/SUBSYSTEM:WINDOWS
)

# Dependencies
target_link_libraries(${PROJECT_NAME} PUBLIC "$<$<CONFIG:Release>:Python3>")
target_link_directories(${PROJECT_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/Python3/libs")
install(TARGETS ${PROJECT_NAME} DESTINATION ${CMAKE_INSTALL_PREFIX})
