cmake_minimum_required(VERSION 3.30.0)
project(mycalc)

add_subdirectory(math)

add_executable(mycalc mycalc.cpp)
target_link_libraries(mycalc SHARED mathlib)