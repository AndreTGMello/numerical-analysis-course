#bin/bash!
if [ -z "$1" ]
  then
    arg="-r"
  else
    arg="$1"
fi
if [ "$arg" = "-r" ] || [ "$arg" = "--run" ];
  then
    #gfortran ./lib/ogpf.f90 num_integration.f95 bspline.f95 utils.f95 g_elimination.f95 main_old.f95 -o main
    #./main
    gfortran ./lib/ogpf.f90 utils.f95 user_func.f95 g_elimination.f95 num_integration.f95 bspline.f95 finite_elements.f95 main.f95 -o main
    ./main
   
  elif [ "$arg" = "-c" ] || [ "$arg" = "--compile" ];
    then
    gfortran ./lib/ogpf.f90 utils.f95 user_func.f95 g_elimination.f95 num_integration.f95 bspline.f95 finite_elements.f95 main.f95 -o main

  elif [ "$arg" = "-t" ] || [ "$arg" = "--test" ] ;
    then
    arg="$2"
    if [ "$arg" = "gauss" ];
      then
        # TESTING GAUSS ELIMINATION
        gfortran g_elimination.f95 utils.f95 ./tests/g_elimination_tests.f95 -o ./tests/g_elimination_tests.test
        ./tests/bm_ge_tests.test

    elif [ "$arg" = "integration" ];
      then
        # TESTING NUMERICAL INTEGRATION
        gfortran utils.f95 num_integration.f95 ./lib/quadpack_double.f90 ./tests/num_integration_tests.f95 -o ./tests/num_integration_tests.test
        ./tests/num_integration_tests.test
  
    elif [ "$arg" = "splines" ];
      then
        # TESTING BSPLINES
        gfortran utils.f95 ./lib/ogpf.f90 user_func.f95 g_elimination.f95 num_integration.f95 finite_elements.f95 bspline.f95 ./tests/bspline_tests.f95 -o ./tests/bspline_tests.test 
        ./tests/bspline_tests.test
  
    elif [ "$arg" = "plot" ];
      then
        # PLOTTING LIBRARY DEMO
        gfortran ./lib/ogpf.f90 ./lib/demo.f90 -o ./lib/plot.test
        ./lib/plot.test
  
    elif [ "$arg" = "valgrind" ];
      then
        # VALGRIND
    	gfortran ./lib/ogpf.f90 utils.f95 user_func.f95 g_elimination.f95 num_integration.f95 bspline.f95 finite_elements.f95 main.f95 -o main
    	gfortran ./lib/ogpf.f90 utils.f95 user_func.f95 g_elimination.f95 num_integration.f95 bspline.f95 finite_elements.f95 main.f95 -g main
        valgrind ./main
    elif [ "$arg" = "gdb" ];
      then
        gfortran ./lib/ogpf.f90 utils.f95 user_func.f95 g_elimination.f95 num_integration.f95 bspline.f95 finite_elements.f95 main.f95 -g -Wall -Wextra -Warray-temporaries -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan main
        gdb ./main
    fi
fi
