# Add test executables
file(GLOB TEST_FILES "*.cpp")

foreach(test_file ${TEST_FILES})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
        "/opt/homebrew/include/eigen3"
    )
    # Link against the source libraries
    target_link_libraries(${test_name} PRIVATE
        ${PROJECT_NAME}  # Main library
        Eigen3::Eigen
        GTest::gtest
        GTest::gtest_main
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
