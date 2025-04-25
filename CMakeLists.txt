cmake_minimum_required(VERSION 3.30)
project(PreciseDotProduct)

set(CMAKE_CXX_STANDARD 20)

find_package(Boost REQUIRED)

add_executable(PreciseDotProduct main.cc)
target_link_libraries(PreciseDotProduct
        PRIVATE
        Boost::boost)
