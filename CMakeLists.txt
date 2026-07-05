cmake_minimum_required(VERSION 3.16)
project(madgwick_orientation_filter VERSION 1.0.0 LANGUAGES CXX)

option(MADGWICK_BUILD_TESTS "Build unit tests" ON)
option(MADGWICK_FETCH_DEPENDENCIES "Download dependencies not found locally" ON)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    if(NOT MADGWICK_FETCH_DEPENDENCIES)
        message(FATAL_ERROR "Eigen3 was not found; install it or enable MADGWICK_FETCH_DEPENDENCIES")
    endif()
    include(FetchContent)
    FetchContent_Declare(
        eigen
        GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
        GIT_TAG 3.4.0
        GIT_SHALLOW TRUE)
    set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
    set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(eigen)
endif()

add_library(madgwick Madgwick.cpp)
add_library(Madgwick::madgwick ALIAS madgwick)
target_include_directories(madgwick PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>)
target_link_libraries(madgwick PUBLIC Eigen3::Eigen)
target_compile_features(madgwick PUBLIC cxx_std_17)

if(MSVC)
    target_compile_options(madgwick PRIVATE /W4 /permissive-)
else()
    target_compile_options(madgwick PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(MADGWICK_BUILD_TESTS)
    set(BUILD_TESTING ON CACHE BOOL "" FORCE)
    include(CTest)
    include(FetchContent)
    FetchContent_Declare(
        googletest
        GIT_REPOSITORY https://github.com/google/googletest.git
        GIT_TAG v1.15.2
        GIT_SHALLOW TRUE)
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_executable(madgwick_tests tests/MadgwickTests.cpp)
    target_link_libraries(madgwick_tests PRIVATE Madgwick::madgwick GTest::gtest_main)
    target_compile_features(madgwick_tests PRIVATE cxx_std_17)
    include(GoogleTest)
    gtest_discover_tests(madgwick_tests)
endif()
