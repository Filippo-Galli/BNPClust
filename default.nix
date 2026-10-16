{
  stdenv,
  cmake,
  eigen,
  R,
  lib,
  rPackages,
  python3Packages,
  zlib,
  zstd,
  xz,
  bzip2,
  libdeflate,
  icu,
  pcre2,
  Rcpp ? rPackages.Rcpp,
  RcppEigen ? rPackages.RcppEigen,
  nanobind ? python3Packages.nanobind,
}:

stdenv.mkDerivation {
  pname = "BNPClust";
  version = "0.0.1-unstable-2026-09-30";

  src = lib.fileset.toSource {
    root = ./.;
    fileset = lib.fileset.unions [
      ./BNPClust
      ./CMakeLists.txt
    ];
  };

  nativeBuildInputs = [
    cmake
    R
    python3Packages.python
  ];

  buildInputs = [
    eigen

    # necessary for python bindings
    nanobind

    # necessary for R bindings
    Rcpp
    RcppEigen

    # Libraries referenced by R
    zlib
    zstd
    xz
    bzip2
    libdeflate
    icu
    pcre2
  ];

  cmakeFlags = [
    "-DBNPCLUST_BUILD_R=ON"
    "-DBNPCLUST_BUILD_PYTHON=ON"

    "-DRCPP_INCLUDE_DIR=${Rcpp}/library/Rcpp/include"
    "-DRCPPEIGEN_INCLUDE_DIR=${RcppEigen}/library/RcppEigen/include"

    "-DPython_EXECUTABLE=${python3Packages.python}/bin/python3"
    "-Dnanobind_DIR=${nanobind}/${python3Packages.python.sitePackages}/nanobind/cmake"

    "-DBNPCLUST_PYTHON_INSTALL_DIR=${python3Packages.python.sitePackages}"
  ];

  meta = {
    description = "Bayesian nonparametric clustering library";
    homepage = "https://github.com/Filippo-Galli/BNPClust";
    license = lib.licenses.gpl3Only;
    platforms = lib.platforms.linux;
  };
}
