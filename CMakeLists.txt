cmake_minimum_required(VERSION 3.16)
project(madgwick_orientation_filter VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(MADGWICK_BUILD_TESTS "Build unit tests" OFF)

add_subdirectory(third_party/eigen EXCLUDE_FROM_ALL)

add_library(madgwick Madgwick.cpp)
target_include_directories(madgwick PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(madgwick PUBLIC Eigen3::Eigen)

if(MADGWICK_BUILD_TESTS)
    include(CTest)
    find_package(GTest REQUIRED)

    add_executable(madgwick_tests tests/MadgwickTests.cpp)
    target_link_libraries(madgwick_tests PRIVATE madgwick GTest::gtest_main)
    include(GoogleTest)
    gtest_discover_tests(madgwick_tests)
endif()
