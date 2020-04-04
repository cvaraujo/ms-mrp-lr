cmake_minimum_required(VERSION 3.14)
project(MaxServiceLagrangean)

set(CMAKE_CXX_STANDARD 11)

include_directories(inc src)

add_executable(MaxServiceLagrangean
        inc/Arc.h
        inc/Graph.h
        inc/Model.h
        src/Arc.cpp
        src/Graph.cpp
        src/Model.cpp
        inc/edmonds_optimum_branching.hpp
        inc/edmonds_optimum_branching_impl.hpp
        main.cpp src/Lagrangean.cpp inc/Lagrangean.h)
