cmake_minimum_required(VERSION 3.14 FATAL_ERROR)


find_package(Eigen3 REQUIRED)

include(FetchContent)
FetchContent_Declare(
        polyscope
        GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
        GIT_TAG        master
)
FetchContent_MakeAvailable(polyscope)


add_executable(generate_cylinder
        executables/generate_cylinder.cc
)

target_include_directories(generate_cylinder PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(generate_cylinder
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        polyscope
)

set_target_properties(generate_cylinder PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)


add_executable(generate_prism
        executables/generate_prism.cc
)

target_include_directories(generate_prism PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(generate_prism
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        polyscope
)

set_target_properties(generate_prism PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)



add_executable(generate_adjacent_prisms
        executables/generate_adjacent_prisms.cc
)

target_include_directories(generate_adjacent_prisms PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(generate_adjacent_prisms
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        polyscope
)

set_target_properties(generate_adjacent_prisms PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)


add_executable(generate_spheritorus
        executables/generate_spheritorus.cc
        header/ovm_nd_types.hh
        header/base_shape_generator.hh
        header/cartesian_product.hh
        header/nd_subdivision.hh
        header/orientation.hh
        validation/validation.hh
        header/homology.hh
        header/projection.hh
        validation/quality_metrics.hh
        validation/visualize_4d.hh
        validation/results_dir.hh
)

target_include_directories(generate_spheritorus PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(generate_spheritorus
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        Eigen3::Eigen
        polyscope
)

target_compile_definitions(generate_spheritorus
        PRIVATE
        RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/results"
)

set_target_properties(generate_spheritorus PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)


add_executable(schwarz_lantern
        validation/schwarz_lantern.cc
        header/ovm_nd_types.hh
        header/base_shape_generator.hh
        header/cartesian_product.hh
        header/nd_subdivision.hh
        validation/validation.hh
        header/projection.hh
        validation/results_dir.hh
)

target_include_directories(schwarz_lantern PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(schwarz_lantern
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
)

target_compile_definitions(schwarz_lantern
        PRIVATE
        RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/results"
)

set_target_properties(schwarz_lantern PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)


add_executable(generate_ditorus
        executables/generate_ditorus.cc
        header/ovm_nd_types.hh
        header/base_shape_generator.hh
        header/cartesian_product.hh
        header/nd_subdivision.hh
        header/orientation.hh
        validation/validation.hh
        header/homology.hh
        header/projection.hh
        validation/quality_metrics.hh
        validation/visualize_4d.hh
        validation/results_dir.hh
)

target_include_directories(generate_ditorus PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(generate_ditorus
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        Eigen3::Eigen
        polyscope
)

target_compile_definitions(generate_ditorus
        PRIVATE
        RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/results"
)

set_target_properties(generate_ditorus PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)


add_executable(mesh_quality_assessment
        validation/mesh_quality_assessment.cc
        header/ovm_nd_types.hh
        header/base_shape_generator.hh
        header/cartesian_product.hh
        header/nd_subdivision.hh
        header/projection.hh
        validation/quality_metrics.hh
        validation/results_dir.hh
)

target_include_directories(mesh_quality_assessment PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(mesh_quality_assessment
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        Eigen3::Eigen
)

target_compile_definitions(mesh_quality_assessment
        PRIVATE
        RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/results"
)

set_target_properties(mesh_quality_assessment PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)


add_executable(timing
        validation/timing.cc
        header/ovm_nd_types.hh
        header/base_shape_generator.hh
        header/cartesian_product.hh
        header/nd_subdivision.hh
        header/orientation.hh
        validation/validation.hh
        header/homology.hh
        header/projection.hh
        validation/results_dir.hh
)

target_include_directories(timing PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/header
        ${CMAKE_CURRENT_SOURCE_DIR}/validation
)

target_link_libraries(timing
        PRIVATE
        OpenVolumeMesh::OpenVolumeMesh
        Eigen3::Eigen
)

target_compile_definitions(timing
        PRIVATE
        RESULTS_DIR="${CMAKE_CURRENT_SOURCE_DIR}/validation/results"
        NDEBUG
)

target_compile_options(timing PRIVATE -O2)

set_target_properties(timing PROPERTIES
        CXX_STANDARD 17
        CXX_STANDARD_REQUIRED YES
        CXX_EXTENSIONS NO
)
