cmake_minimum_required(VERSION 3.13)
project(Comparison)

INCLUDE_DIRECTORIES($ENV{HOME} $ENV{HOME}/libscapi/install/include)
link_directories(/usr/ssl/lib/ $ENV{HOME}/libscapi/install/lib)

set(CMAKE_CXX_STANDARD 14)

add_executable(Comparison main.cpp ProtocolParty.cpp ProtocolParty.h)


TARGET_LINK_LIBRARIES(Comparison ntl.a
        $ENV{HOME}/libscapi/libscapi.a gmp
        boost_system pthread crypto dl)
