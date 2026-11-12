cmake_minimum_required(VERSION 3.16)
project(gpxsmooth VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_AUTOMOC ON)

find_package(Qt6 REQUIRED COMPONENTS Widgets)

add_executable(gpxsmooth
    src/main.cpp
    src/gpxtrack.cpp
    src/gpxtrack.h
    src/graphwidget.cpp
    src/graphwidget.h
)

target_link_libraries(gpxsmooth PRIVATE Qt6::Widgets)
