cmake_minimum_required(VERSION 3.5)
project(handeyecat)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
add_compile_options(-Wall -Wextra -Wpedantic)

find_package(Eigen3 REQUIRED)

add_library(${PROJECT_NAME} src/handeyecat.cpp)
target_include_directories(${PROJECT_NAME} SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIRS})

add_executable(${PROJECT_NAME}_demo src/demo.cpp)
target_link_libraries(${PROJECT_NAME}_demo ${PROJECT_NAME})
