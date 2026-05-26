cmake_minimum_required(VERSION 3.16)
project(BNPClust LANGUAGES CXX)

# --- Standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# --- Default to Release if nothing specified ---
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# --- Sources ---
file(GLOB_RECURSE BNPCLUST_SOURCES CONFIGURE_DEPENDS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")
list(FILTER BNPCLUST_SOURCES EXCLUDE REGEX ".*/r_bindings\\.cpp$")

add_library(bnpclust_core SHARED ${BNPCLUST_SOURCES})
target_include_directories(bnpclust_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/src")

# --- Eigen ---
find_package(Eigen3 REQUIRED)
target_link_libraries(bnpclust_core PUBLIC Eigen3::Eigen)
target_compile_definitions(bnpclust_core PRIVATE
    EIGEN_NO_DEBUG
    EIGEN_DONT_PARALLELIZE   # let OpenMP handle parallelism instead
)

# --- OpenMP ---
find_package(OpenMP REQUIRED)
target_link_libraries(bnpclust_core PUBLIC OpenMP::OpenMP_CXX)

# --- Machine-specific flags (only when building natively, not cross-compiling) ---
if(NOT CMAKE_CROSSCOMPILING)
  target_compile_options(bnpclust_core PRIVATE
    $<$<CONFIG:Release>:-march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize>
    $<$<CONFIG:Debug>:-fsanitize=address -fno-omit-frame-pointer>
  )
  target_link_options(bnpclust_core PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
  )
endif()

set_target_properties(bnpclust_core PROPERTIES OUTPUT_NAME "bnpclust_core")
