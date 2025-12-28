cmake_minimum_required(VERSION 3.31)

project(raytrace LANGUAGES CXX)

add_executable(raytrace
        vec3.hpp
        color.hpp
        ray.hpp)
target_sources(raytrace PUBLIC main.cpp)
