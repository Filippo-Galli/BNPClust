{
  mkShell,
  rWrapper,
  rPackages,
  doxygen,
  aria2,
  mcclustExt,
}:

mkShell {
  packages = [
    doxygen
    aria2
    # Install R with the packages
    (rWrapper.override {
      packages = with rPackages; [
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

      ];
    })
  ];
}
