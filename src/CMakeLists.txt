# Add individual library targets for each module first
add_library(Laminate STATIC
    Laminate/cltweight.cpp
    Laminate/laminateSection.cpp
    Laminate/lpfeasible.cpp
    Laminate/varsize.cpp
)

add_library(Section STATIC
    Section/section.cpp
    Section/sectionBase.cpp
)

add_library(GlobalOptimiser STATIC
    GlobalOptimiser/approxFunction.cpp
    GlobalOptimiser/globalOptimiser.cpp
    GlobalOptimiser/scminmaxProb.cpp
)

add_library(InterfaceOptimiser STATIC
    InterfaceOptimiser/dampingInterface.cpp
    InterfaceOptimiser/responseInterface.cpp
)

add_library(SDPA STATIC
    SDPA/parameter.cpp
    SDPA/scalarSDP.cpp
    SDPA/SDPA.cpp
)

add_library(BoundSDP STATIC
    BoundSDP/boundSDP.cpp
)

add_library(TestProblems STATIC
    TestProblems/designProblem1.cpp
    TestProblems/solutionProblem1.cpp
)

# Add source files for the main library
file(GLOB_RECURSE SOURCES 
    "BoundSDP/*.cpp"
    "GlobalOptimiser/*.cpp"
    "InterfaceOptimiser/*.cpp"
    "Laminate/*.cpp"
    "Miki/*.cpp"
    "SDPA/*.cpp"
    "Section/*.cpp"
    "TestProblems/*.cpp"
    "fSDPCone/*.cpp"
)

# Create the main library
add_library(${PROJECT_NAME} STATIC ${SOURCES})

# Include directories for main library
target_include_directories(${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    "/opt/homebrew/include/eigen3"
)

# Link the main library with its dependencies
target_link_libraries(${PROJECT_NAME} PUBLIC
    Laminate
    Section
    GlobalOptimiser
    InterfaceOptimiser
    SDPA
    BoundSDP
    TestProblems
    Eigen3::Eigen
)

# Make headers available for the tests
target_include_directories(Laminate PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(Section PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(GlobalOptimiser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(InterfaceOptimiser PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(SDPA PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(BoundSDP PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_include_directories(TestProblems PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
