cmake_minimum_required(VERSION 3.30.0)
project(mathlib)

file(GLOB SOURCES "*.cpp")
add_library(mathlib SHARED ${SOURCES})