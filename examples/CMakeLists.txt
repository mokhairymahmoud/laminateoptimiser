set(ACTIVE_EXAMPLE_FILES
    mock_calculix_e2e.cpp
)

foreach(example_file ${ACTIVE_EXAMPLE_FILES})
    get_filename_component(example_name ${example_file} NAME_WE)
    add_executable(${example_name} ${example_file})
    target_include_directories(${example_name} PRIVATE
        ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(${example_name} PRIVATE
        ${PROJECT_NAME}
        Eigen3::Eigen
    )
endforeach()
