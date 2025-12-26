cmake_minimum_required(VERSION 4.0)
project(src)

set(CMAKE_CXX_STANDARD 11)

add_executable(src main.cpp
        matrix.h
        matrix.cpp
)
