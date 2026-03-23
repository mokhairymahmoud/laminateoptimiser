set(ACTIVE_TEST_FILES
    Miki_test.cpp
    calculix_backend_test.cpp
    dampingUtils_test.cpp
    global_driver_test.cpp
    lamination_parameter_derivatives_test.cpp
    laminate_test.cpp
    minmaxprob_test.cpp
    section_test.cpp
    thickness_test.cpp
)

foreach(test_file ${ACTIVE_TEST_FILES})
    get_filename_component(test_name ${test_file} NAME_WE)
    add_executable(${test_name} ${test_file})
    target_include_directories(${test_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(${test_name} PRIVATE
        ${PROJECT_NAME}
        Eigen3::Eigen
        GTest::gtest
        GTest::gtest_main
    )
    add_test(NAME ${test_name} COMMAND ${test_name})
endforeach()
