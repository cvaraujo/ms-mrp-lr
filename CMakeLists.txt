cmake_minimum_required(VERSION 3.5)
project(MSLagrangean)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DIL_STD -D_GLIBCXX_USE_CXX11_ABI=0 -std=c++0x ")

include_directories(headers src /opt/gurobi902/linux64/include /opt/boost_1_73_0)

link_directories(/opt/gurobi902/linux64/lib /opt/boost_1_73_0)

set(SOURCE_FILES
        headers/Arc.h
        headers/Graph.h
        headers/Model.h
        headers/BarrierMethod.h
        src/Arc.cpp
        src/Graph.cpp
        src/Model.cpp
        src/BarrierMethod.cpp
        headers/edmonds_optimum_branching.hpp
        headers/edmonds_optimum_branching_impl.hpp
        main.cpp src/Lagrangean.cpp headers/Lagrangean.h)

add_executable(MSLagrangean ${SOURCE_FILES})

target_link_libraries(MSLagrangean gurobi_c++ gurobi90 pthread m)
