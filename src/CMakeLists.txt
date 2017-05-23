project(triangles)
cmake_minimum_required(VERSION 2.6)

set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_SOURCE_DIR}/..)
set(INSTALL_NAME_DIR ${CMAKE_INSTALL_PREFIX}/bin)
set(SRC_DIR ${PROJECT_SOURCE_DIR})

list(APPEND SRC ${SRC_DIR}/main.cpp)

add_definitions(-Wall -Werror -std=c++14)

add_executable(${PROJECT_NAME} ${SRC})

target_link_libraries(${PROJECT_NAME})

install(TARGETS ${PROJECT_NAME} DESTINATION ${INSTALL_NAME_DIR})