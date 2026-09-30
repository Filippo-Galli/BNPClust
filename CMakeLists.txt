cmake_minimum_required(VERSION 3.16)
project(BNPClust LANGUAGES CXX)

# --- Standard ---
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_compile_options(-fno-finite-math-only)

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
    $<$<CONFIG:Release>:-march=native -mtune=native -O3 -funroll-loops -ftree-vectorize>
    $<$<CONFIG:Debug>:-fsanitize=address -fno-omit-frame-pointer>
  )
  target_link_options(bnpclust_core PRIVATE
    $<$<CONFIG:Debug>:-fsanitize=address>
  )
endif()

set_target_properties(bnpclust_core PROPERTIES OUTPUT_NAME "bnpclust_core")

##################################
# R / Rcpp / RcppEigen
##################################
find_package(PkgConfig QUIET)
if(PkgConfig_FOUND)
  pkg_check_modules(R IMPORTED_TARGET libR)
endif()

if(NOT TARGET PkgConfig::R)
  # Fallback to R CMD config if pkg-config is not available
  find_program(R_EXECUTABLE R REQUIRED)
  execute_process(
    COMMAND ${R_EXECUTABLE} CMD config --cppflags
    OUTPUT_VARIABLE R_CPPFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  execute_process(
    COMMAND ${R_EXECUTABLE} CMD config --ldflags
    OUTPUT_VARIABLE R_LDFLAGS
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  separate_arguments(R_CPPFLAGS)
  separate_arguments(R_LDFLAGS)
endif()

# Allow passing RCPP_INCLUDE_DIR and RCPPEIGEN_INCLUDE_DIR from CMake command line (-D...)
# This is crucial for Nix derivations to avoid running Rscript in pure build sandboxes
if(NOT RCPP_INCLUDE_DIR)
  find_program(RSCRIPT_EXECUTABLE Rscript REQUIRED)
  execute_process(
    COMMAND env -u R_LIBS_USER ${RSCRIPT_EXECUTABLE} -e "cat(system.file('include', package='Rcpp'))"
    OUTPUT_VARIABLE RCPP_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

if(NOT RCPPEIGEN_INCLUDE_DIR)
  find_program(RSCRIPT_EXECUTABLE Rscript REQUIRED)
  execute_process(
    COMMAND env -u R_LIBS_USER ${RSCRIPT_EXECUTABLE} -e "cat(system.file('include', package='RcppEigen'))"
    OUTPUT_VARIABLE RCPPEIGEN_INCLUDE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )
endif()

add_library(bnpclust_r SHARED
  src/r_bindings.cpp
)

target_include_directories(bnpclust_r PRIVATE
  "${CMAKE_CURRENT_SOURCE_DIR}/src"
  "${RCPP_INCLUDE_DIR}"
  "${RCPPEIGEN_INCLUDE_DIR}"
)

# ENABLE_LEGACY_NONAPI is needed for R >= 4.6 compatibility when using Rcpp headers that reference R_NamespaceRegistry
target_compile_definitions(bnpclust_r PRIVATE
  ENABLE_LEGACY_NONAPI
)

target_link_libraries(bnpclust_r PRIVATE
  bnpclust_core
)

if(TARGET PkgConfig::R)
  target_link_libraries(bnpclust_r PRIVATE PkgConfig::R)
else()
  target_compile_options(bnpclust_r PRIVATE ${R_CPPFLAGS})
  target_link_options(bnpclust_r PRIVATE ${R_LDFLAGS})
endif()

include(GNUInstallDirs)

install(
  TARGETS bnpclust_core
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

if(BNPCLUST_BUILD_R)
  install(
    TARGETS bnpclust_r
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()

set_target_properties(bnpclust_r PROPERTIES
  OUTPUT_NAME "bnpclust_r"
  BUILD_RPATH "$ORIGIN"
  INSTALL_RPATH "$ORIGIN"
)
