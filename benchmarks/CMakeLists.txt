set(ACTIVE_BENCHMARK_FILES
    calculix_e2e_test.cpp
    calculix_composipy_benchmark_test.cpp
    calculix_laminate_parameter_benchmark_test.cpp
)

foreach(benchmark_file ${ACTIVE_BENCHMARK_FILES})
    get_filename_component(benchmark_name ${benchmark_file} NAME_WE)
    add_executable(${benchmark_name} ${benchmark_file})
    target_include_directories(${benchmark_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(${benchmark_name} PRIVATE
        ${PROJECT_NAME}
        Eigen3::Eigen
        GTest::gtest
        GTest::gtest_main
    )
    add_test(NAME ${benchmark_name} COMMAND ${benchmark_name})
endforeach()
