## Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'RelWithDebInfo' as none was specified.")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_C_STANDARD 11)
# to enable cuda : uncomment
# enable_language(CUDA)

add_compile_options(
    -Wall 
    -Wextra 
    -Wshadow 
    -Wno-unknown-pragmas 
    -Wno-unused-parameter
    -O3
)

file(GLOB APP_SOURCES *.cpp *.c)
message("sources = ${APP_SOURCES}")
add_library(
    femsolver
    ${APP_SOURCES}
)
find_library(m m)
target_link_libraries(femsolver PUBLIC m)
