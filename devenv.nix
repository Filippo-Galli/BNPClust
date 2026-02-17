{
  pkgs,
  lib,
  config,
  inputs,
  ...
}:

let
  openblas = pkgs.openblas.overrideAttrs (old: old // { isILP64 = false; });
  RWithOpenBlas = pkgs.R.override {
    blas = openblas;
    lapack = openblas;
  };
in
{
  env.GREET = "Good day, happy coding";

  # ========== DISABLE NIX NATIVE ENFORCEMENT FOR PERFORMANCE ==========
  env.NIX_ENFORCE_NO_NATIVE = "0";

  # ========== OPTIMIZED C++ COMPILATION FLAGS FOR MCMC ==========
  # CPU architecture-specific optimizations (best performance)
  env.MARCH_FLAGS = "-march=native -mtune=native";

  # Optimization levels for MCMC (computational intensity is high)
  # -O3: aggressive optimizations
  # -ffast-math: allow aggressive floating point optimizations (safe for MCMC)
  # -funroll-loops: unroll loops for better CPU cache utilization
  # -ftree-vectorize: auto-vectorize loops (SIMD)
  # -flto: link-time optimization
  env.CXX_STD = "CXX17";
  env.PKG_CXXFLAGS = "-O2 -g -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize -flto=auto -fopenmp";
  env.CXXFLAGS = "-O2 -g -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize -flto=auto -fopenmp";
  env.PKG_CFLAGS = "-O2 -g -march=native -mtune=native -ffast-math -funroll-loops -ftree-vectorize -flto=auto -fopenmp";

  # Eigen specific flags
  env.PKG_CXXFLAGS_EIGEN = "-DEIGEN_NO_DEBUG -DEIGEN_DONT_PARALLELIZE";

  # Linker flags for optimization and OpenMP
  env.PKG_LIBS = "-L${openblas}/lib -lopenblas -flto=auto -fopenmp";
  env.LDFLAGS = "-L${openblas}/lib -lopenblas -flto=auto -fopenmp";

  # Set R environment variables
  env.R_LIBS_USER = "${config.env.DEVENV_STATE}/R";
  env.R_LIBS_SITE = "${pkgs.R}/library";
  env.PKG_CONFIG_PATH = "${pkgs.pkg-config}/lib/pkgconfig:${pkgs.openssl.dev}/lib/pkgconfig:${pkgs.udunits}/lib/pkgconfig";
  env.LD_LIBRARY_PATH = "${openblas}/lib:${pkgs.openssl}/lib:${pkgs.udunits}/lib:${pkgs.geos}/lib:${pkgs.gdal}/lib:${pkgs.proj}/lib:${pkgs.cairo}/lib:${pkgs.pango}/lib:${pkgs.fontconfig}/lib:${pkgs.libpng}/lib:${pkgs.xorg.libX11}/lib";

  env.LIBRARY_PATH = "${openblas}/lib:${pkgs.openssl.dev}/lib:${pkgs.udunits}/lib:${pkgs.geos}/lib:${pkgs.gdal}/lib:${pkgs.proj}/lib:${pkgs.cairo}/lib:${pkgs.pango}/lib:${pkgs.fontconfig}/lib:${pkgs.libpng}/lib:${pkgs.xorg.libX11}/lib";

  env.OPENBLAS_NUM_THREADS = "1";

  # Podman
  env.PODMAN_USERNS = "keep-id";

  packages = [
    pkgs.git
    pkgs.podman
    pkgs.podman-compose
    pkgs.aria2
    
    # R development with C++ support
    RWithOpenBlas
    pkgs.rPackages.Rcpp
    pkgs.rPackages.RcppEigen

    # Graphics rendering packages for R PNG output
    pkgs.cairo
    pkgs.pango
    pkgs.fontconfig
    pkgs.xorg.libX11
    pkgs.libpng

    # VS Code R integration packages
    # pkgs.rPackages.httpgd # download it manually due to broken nixpkgs version
    pkgs.rPackages.languageserver
    pkgs.rPackages.jsonlite
    pkgs.rPackages.renv
    pkgs.rPackages.mvtnorm
    pkgs.rPackages.gtools

    # Data analysis packages
    pkgs.rPackages.ggplot2
    pkgs.rPackages.dplyr
    pkgs.rPackages.tidyr

    # MCMC and clustering tools
    pkgs.rPackages.spam
    pkgs.rPackages.fields
    pkgs.rPackages.viridisLite
    pkgs.rPackages.RColorBrewer
    pkgs.rPackages.pheatmap
    pkgs.rPackages.mcclust
    pkgs.rPackages.salso
    pkgs.rPackages.mclust

    # ========== ADD THESE FOR SPATIAL PACKAGES (s2, sf, units, spdep, rgeoda) ==========

    # System dependencies for spatial packages
    pkgs.openssl # Required by s2
    pkgs.udunits # Required by units
    pkgs.geos # Required by sf (s2 internally)
    pkgs.gdal # Required by sf
    pkgs.proj # Required by sf
    pkgs.libxml2 # Often needed by spatial packages
    pkgs.boost.dev # Required by rgeoda

    # Spatial packages
    pkgs.rPackages.s2
    pkgs.rPackages.sf
    pkgs.rPackages.units
    pkgs.rPackages.spdep
    pkgs.rPackages.rgeoda

    # ========== END SPATIAL PACKAGES ==========

    # C++ development tools
    pkgs.gcc13
    pkgs.gfortran
    pkgs.pkg-config
    pkgs.cmake
    pkgs.llvmPackages_18.clang-tools
    pkgs.perf

    # OpenMP support
    pkgs.llvmPackages_18.openmp

    # Linear algebra libraries (for RcppEigen)
    pkgs.eigen
    openblas

    # Additional R development packages
    pkgs.rPackages.devtools
    pkgs.rPackages.testthat
    pkgs.rPackages.roxygen2
    pkgs.rPackages.knitr
    pkgs.rPackages.aricode
    pkgs.rPackages.cluster
    pkgs.rPackages.reshape2
    pkgs.rPackages.label_switching
    pkgs.rPackages.survival
    pkgs.rPackages.glmnet

    # Documentation
    pkgs.doxygen
  ];

  scripts.r-setup.exec = ''
    echo "Setting up R environment..."
    mkdir -p $R_LIBS_USER
    echo "R library path: $R_LIBS_USER"
    echo "R site library: $R_LIBS_SITE"
    echo "PKG_CONFIG_PATH: $PKG_CONFIG_PATH"
    echo "LD_LIBRARY_PATH: $LD_LIBRARY_PATH"
    echo "You can now use R with Rcpp, RcppEigen, spatial packages, and SALSO support!"
  '';

  scripts.install-mcclust-ext.exec = ''
    echo "Installing mcclust.ext from Warwick archive..."
    mkdir -p $R_LIBS_USER
    chmod -R u+w $R_LIBS_USER

    export R_LIBS=$R_LIBS_USER

    R -e "
      .libPaths(Sys.getenv('R_LIBS_USER'))
      if (!require('devtools', quietly = TRUE)) {
        install.packages('devtools', repos='https://cloud.r-project.org', lib=.libPaths()[1])
      }
      devtools::install_url('http://wrap.warwick.ac.uk/71934/1/mcclust.ext_1.0.tar.gz', lib=.libPaths()[1])
    "
    echo "✅ mcclust.ext installed successfully!"
  '';

  scripts.setup-clangd.exec = ''
        echo "Setting up clangd configuration..."
        
        R_INCLUDE=$(R --slave -e "cat(R.home('include'))")
        RCPP_INCLUDE=$(R --slave -e "cat(system.file('include', package='Rcpp'))")
        RCPPEIGEN_INCLUDE=$(R --slave -e "cat(system.file('include', package='RcppEigen'))")
        
        # Get GCC paths for C++ standard library
        LIBCXX_INCLUDE=${pkgs.gcc13.cc}/include/c++/${pkgs.gcc13.cc.version}
        
        cat > .clangd <<EOF
    CompileFlags:
      Add: [
        "-std=c++17",
        "-I./include",
        "-I${pkgs.eigen}/include/eigen3",
        "-I$R_INCLUDE",
        "-I$RCPP_INCLUDE", 
        "-I$RCPPEIGEN_INCLUDE",
        "-I$LIBCXX_INCLUDE",
        "-I$LIBCXX_INCLUDE/x86_64-unknown-linux-gnu"
      ]
    EOF

        echo "✅ Generated .clangd with dynamic include paths:"
        cat .clangd
  '';

  scripts.clean_o_files.exec = ''
    echo "Cleaning up .o files..."
    find . -name "*.o" -type f -delete
    echo "✅ Removed all .o files"
  '';

  enterShell = ''
        setup-clangd

        # Create containers config directory
        mkdir -p ~/.config/containers
        
        # Create a simple permissive policy.json
        cat > ~/.config/containers/policy.json << 'EOF'
    {
      "default": [
        {
          "type": "insecureAcceptAnything"
        }
      ]
    }
    EOF

        echo ""
        echo "🚀 R + C++ + Spatial development environment is ready for VS Code!"
        echo "   - httpgd: Plot viewing"
        echo "   - languageserver: IntelliSense and code completion"
        echo "   - salso: Modern clustering analysis"
        echo "   - sf/s2/units: Spatial data analysis"
        echo "   - spdep: Spatial dependence analysis"
        echo "   - clangd: C++ LSP with correct R/Rcpp paths"
        echo "   - OpenMP: Parallel computing support enabled"
        echo ""
        echo "💡 Run 'install-mcclust-ext' to install mcclust.ext"
        echo "💡 Run 'clean_o_files' to delete all .o files <=> recompile all"
        echo "💡 Run 'test-spatial' to verify spatial package installations"
        echo "💡 Run 'test-openmp' to verify OpenMP is working"
        echo "💡 Run 'check-r-packages' to check all R packages"
        echo "💡 Run 'r-setup' to display environment paths"
        echo ""
  '';

  languages.r.enable = true;
}
