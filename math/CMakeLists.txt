cmake_minimum_required(VERSION 3.30.0)
project(mathlib)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

file(GLOB SOURCES "*.cpp")
add_library(mathlib STATIC ${SOURCES})