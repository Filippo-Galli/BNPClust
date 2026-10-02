{
  mkShell,
  rWrapper,
  rPackages,
  doxygen,
  aria2,
  mcclustExt,
  cmake,
  which,
}:
let
  s2-fixed = rPackages.s2.overrideAttrs (old: {
    nativeBuildInputs = (old.nativeBuildInputs or [ ]) ++ [
      cmake
      which
    ];
    preConfigure = (old.preConfigure or "") + ''
      export S2_FORCE_BUNDLED_ABSEIL=1
    '';
  });
  rPackages_s2Fixed = rPackages.override {
    overrides = {
      s2 = s2-fixed;
    };
  };
in
mkShell {
  packages = [
    doxygen
    aria2

    # Install R with the packages
    (rWrapper.override {
      packages = with rPackages_s2Fixed; [
        Rcpp
        RcppEigen
        ggplot2
        dplyr
        tidyr
        spam
        fields
        pheatmap
        mcclust
        mcclustExt
        mvtnorm
        gtools
        salso
        aricode
        reshape2
        label_switching
        coda
        s2
        sf
        sn
        argparser
        spdep
        dplyr
        readxl
      ];
    })
  ];
}
