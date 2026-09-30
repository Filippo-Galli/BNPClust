{
  stdenv,
  src,
  cmake,
  eigen,
  R,
  lib,
  rPackages,
  zlib,
  zstd,
  xz,
  bzip2,
  libdeflate,
  icu,
  pcre2,
}:

stdenv.mkDerivation {
  pname = "BNPClust";
  version = "0.0.1-unstable-2026-09-30";

  inherit src;

  nativeBuildInputs = [
    cmake
    R
  ];

  buildInputs = [
    eigen
    rPackages.Rcpp
    rPackages.RcppEigen

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
  ];

  meta = {
    description = "Bayesian nonparametric clustering library";
    homepage = "https://github.com/Filippo-Galli/BNPClust";
    license = lib.licenses.gpl3Only;
    platforms = lib.platforms.linux;
  };
}
