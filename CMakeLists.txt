cmake_minimum_required(VERSION 3.5)
project(MSLagrangean)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(/opt/gurobi903/linux64/include /opt/boost_1_73_0 src/ headers/)

link_directories(/opt/gurobi903/linux64/lib /opt/boost_1_73_0)

set(SOURCE_FILES
        main.cpp
        src/Graph.cpp
        headers/Graph.h
        src/Model.cpp
        headers/Model.h
        src/Arc.cpp
        headers/Arc.h
        headers/BarrierMethod.h
        headers/Lagrangean.h
        src/Lagrangean.cpp
        src/BarrierMethod.cpp
        headers/edmonds_optimum_branching.hpp
        headers/edmonds_optimum_branching_impl.hpp)

add_executable(MSLagrangean ${SOURCE_FILES})

target_link_libraries(MSLagrangean gurobi_c++ gurobi90 pthread m)
