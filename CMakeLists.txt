cmake_minimum_required(VERSION 3.5)
project(gec)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -O2")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -O2")

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")


add_subdirectory(irbem)

add_library(cnpy SHARED "cnpy/cnpy.cpp")

add_executable(gec
        main.cpp
        std_atm.h
        geomagnetic.cpp
        geomagnetic.h
        conductivity.h
        )

enable_language(Fortran)

target_link_libraries(gec
        PUBLIC ${Boost_LIBRARIES}
        irbem
        cnpy
        z
        )

target_include_directories(gec SYSTEM PUBLIC ${Boost_INCLUDE_DIR})
