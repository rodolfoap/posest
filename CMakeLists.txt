CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(POSEST)

INCLUDE_DIRECTORIES(mlsl)

ADD_SUBDIRECTORY(levmar)

ADD_LIBRARY(posest buckets.c lqs.c ransac.c rngs.c prosac.c deal.c p3p.c p4pf.c planep4p.c svd3.c polysolve.c
  poseproj.c posest.c align.c lhm.c sam.c compiler.h lqs.h ransac.h prosac.h p3p.h p4pf.h planep4p.h svd3.h
  polysolve.h poseproj.h posest.h util.h rngs.h sam.h
  mlsl/mlsl.c mlsl/mt19937ar.c mlsl/redblack.c mlsl/sobolseq.c
)
TARGET_INCLUDE_DIRECTORIES(posest PRIVATE levmar)

OPTION(BUILD_DEMOS "Build demo programs?" TRUE)

# demo program
IF(BUILD_DEMOS)
  # paths to levmar & lapack/blas libraries

  SET(LIBS posest levmar m)
  SET(LIBS ${LIBS} ${LAPACKBLAS_LIB_NAMES} ${F2C_LIB_NAME})

  #LINK_DIRECTORIES(${CMAKE_BINARY_DIR}) # location of the posest library
  #LINK_DIRECTORIES(${LAPACKBLAS_DIR})

  ADD_EXECUTABLE(posest_demo posest_demo.c)
  TARGET_LINK_LIBRARIES(posest_demo ${LIBS})

  ADD_EXECUTABLE(binocposest_demo binocposest_demo.c)
  TARGET_LINK_LIBRARIES(binocposest_demo ${LIBS})

  MESSAGE(STATUS "demos will be linked against ${LIBS}")

# make sure that the libraries are built before the demos
  ADD_DEPENDENCIES(posest_demo posest levmar)
  ADD_DEPENDENCIES(binocposest_demo posest levmar)
ENDIF(BUILD_DEMOS)
