# tests/CMakeLists.txt
cmake_minimum_required(VERSION 3.10)
project(BamToolsTests)

# If you haven't already in the top-level CMakeLists:
enable_testing()

# Find Google Test (assuming it's installed system-wide and provides a GTestConfig.cmake or gtest.pc)
find_package(GTest REQUIRED)
find_package(Threads REQUIRED)

# Include Google Test's CMake functions (this might be optional in some installs)
include(GoogleTest)

# Add test target
add_executable(testall 
   main_test.cpp
   test_example.cpp
   test_cut.cpp
   test_valid.cpp
)

# Link GTest libraries and your BamTools library
target_link_libraries(testall PRIVATE 
   GTest::gtest_main 
   BamTools
   Threads::Threads
)

# Discover and register tests with CTest
gtest_discover_tests(testall)

