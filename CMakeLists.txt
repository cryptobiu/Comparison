cmake_minimum_required(VERSION 3.10)
project(Comparison)

INCLUDE_DIRECTORIES($ENV{HOME} $ENV{HOME}/libscapi/install/include)
link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -mbmi2")

add_executable(Comparison main.cpp ProtocolParty.cpp ProtocolParty.h SmallFields.cpp SmallFields.h)


TARGET_LINK_LIBRARIES(Comparison ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp
        boost_system pthread crypto dl)
