set(ACTIVE_SOURCES
    BoundSDP/boundSDP.cpp
    GlobalOptimiser/approxFunction.cpp
    GlobalOptimiser/scminmaxProb.cpp
    Laminate/cltweight.cpp
    Laminate/laminateSection.cpp
    Laminate/lpfeasible.cpp
    Laminate/varsize.cpp
    Miki/Miki.cpp
    SDPA/parameter.cpp
    SDPA/scalarSDP.cpp
    SDPA/SDPA.cpp
    Section/section.cpp
    Section/sectionBase.cpp
)

add_library(${PROJECT_NAME} STATIC ${ACTIVE_SOURCES})

target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(${PROJECT_NAME} PUBLIC
    Eigen3::Eigen
)
