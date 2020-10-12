# Invoking the CMake build from the command line is a two step process, first, *generate* a build by running the following:
#
#   cmake -Bbuild -H. -DCMAKE_BUILD_TYPE=Release
# 
# Where `Release` can be replaced with `Debug`, which contains gdb and address sanitiser definitions already written for you.
# Once a build is created, proceed with *compilation*:
# 
#   cmake --build build --target cg_template --config Release # optionally, for parallel build, append -j $(nproc)
#
# This creates the executable in the build directory. You only need to *generate* a build if you modify the CMakeList.txt file.
# For any other changes to the source code, simply recompile.


cmake_minimum_required(VERSION 3.12)
project(RedNoise)

set(CMAKE_CXX_STANDARD 14)

list(APPEND CMAKE_PREFIX_PATH libs/glm-0.9.7.2/cmake)

find_package(glm REQUIRED)
find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS} ${GLM_INCLUDE_DIRS})

add_executable(RedNoise
        libs/sdw/CanvasPoint.cpp
        libs/sdw/CanvasTriangle.cpp
        libs/sdw/Colour.cpp
        libs/sdw/DrawingWindow.cpp
        libs/sdw/ModelTriangle.cpp
        libs/sdw/RayTriangleIntersection.cpp
        libs/sdw/TextureMap.cpp
        libs/sdw/TexturePoint.cpp
        libs/sdw/Utils.cpp
        src/RedNoise.cpp)

target_compile_options(RedNoise
        PUBLIC
        -Wall
        -Wextra
        -Wcast-align
        -Wfatal-errors
        -Werror=return-type
        -Wno-unused-parameter
        -Wno-unused-variable
        -Wno-ignored-attributes)

set(DEBUG_OPTIONS -O2 -fno-omit-frame-pointer -g)
set(RELEASE_OPTIONS -O3 -march=native -mtune=native)

target_link_libraries(RedNoise PUBLIC $<$<CONFIG:Debug>:-Wl,-lasan>)

include_directories(libs/sdw)

target_compile_options(RedNoise PUBLIC "$<$<CONFIG:RelWithDebInfo>:${RELEASE_OPTIONS}>")
target_compile_options(RedNoise PUBLIC "$<$<CONFIG:Release>:${RELEASE_OPTIONS}>")
target_compile_options(RedNoise PUBLIC "$<$<CONFIG:Debug>:${DEBUG_OPTIONS}>")

target_link_libraries(RedNoise PUBLIC ${GLM_LIBRARIES} ${SDL2_LIBRARIES})
