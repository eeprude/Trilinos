INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../unit_tests)

# first ensure that these cmake boolean variables are defined
ASSERT_DEFINED(
  ${PACKAGE_NAME}_ENABLE_Epetra
  ${PACKAGE_NAME}_ENABLE_Galeri
  ${PACKAGE_NAME}_ENABLE_Isorropia
  ${PACKAGE_NAME}_ENABLE_Tpetra
  )
  
set(XML_FILES_ALL
  challenges/Elasticity2D_reference.xml
  challenges/Elasticity2D_parameters.xml
  challenges/Recirc2D_reference.xml
  challenges/Recirc2D_parameters.xml
  blocked_rebalancing.xml
  ml_ParameterList.xml
  ml_RebalancedParameterList.xml
  muelu_ParameterList.xml
  n1_easy.xml
  n1_easy_pg.xml
  n2_easy_agg.xml
  n2_easy_agg2.xml
  n2_easy_agg3.xml
  n2_easy_export.xml
  s1_easy.xml
  s1_easy_3levels_smoothed.xml
  s1_easy_3levels_unsmoothed.xml
  s1_easy_exercise.xml
  s1_easy_jacobi.xml
  s1_easy_jacobi2.xml
  s2a.xml
  s2_adv_a.xml
  s2_adv_b.xml
  s2_adv_c.xml
  s2b.xml
  s2_easy.xml
  s3a.xml
  s3a_easy.xml
  s3b.xml
  s3b1.xml
  s3b2.xml
  s3b3.xml
  s3b_easy.xml
  s3c.xml
  s3c_easy.xml
  s4a.xml
  s4av.xml
  s4b.xml
  s4c.xml
  s5a.xml
  s5b.xml
  s6_dep.xml
  s6_export.xml
)

TRIBITS_COPY_FILES_TO_BINARY_DIR(CopyXmlFiles
  SOURCE_FILES
    blocked_rebalancing.xml
    ml_ParameterList.xml
    ml_RebalancedParameterList.xml
    muelu_ParameterList.xml
    n1_easy.xml
    n1_easy_pg.xml
    n2_easy_agg.xml
    n2_easy_agg2.xml
    n2_easy_agg3.xml
    n2_easy_export.xml
    s1_easy.xml
    s1_easy_3levels_smoothed.xml
    s1_easy_3levels_unsmoothed.xml
    s1_easy_exercise.xml
    s1_easy_jacobi.xml
    s1_easy_jacobi2.xml
    s2a.xml
    s2_adv_a.xml
    s2_adv_b.xml
    s2_adv_c.xml
    s2b.xml
    s2_easy.xml
    s3a.xml
    s3a_easy.xml
    s3b.xml
    s3b1.xml
    s3b2.xml
    s3b3.xml
    s3b_easy.xml
    s3c.xml
    s3c_easy.xml
    s4a.xml
    s4av.xml
    s4b.xml
    s4c.xml
    s5a.xml
    s5b.xml
    s6_dep.xml
    s6_export.xml
    hands-on.py
  )

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_MLParameterList
    SOURCES MLParameterList.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Tutorial_MLParameterList
    NAME "Test_Tutorial_MLParameterList"
    ARGS "--xml=ml_ParameterList.xml"
    COMM serial mpi
    NUM_MPI_PROCS 4
    )

  # TRIBITS_ADD_TEST(
  #   Tutorial_MLParameterList
  #   NAME "Test_Tutorial_MLParameterList_Rebalanced"
  #   ARGS "--xml=ml_RebalancedParameterList.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 4
  #   )

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_ScalingTest
    SOURCES ScalingTest.cpp
    COMM serial mpi
    )

  TRIBITS_ADD_TEST(
    Tutorial_ScalingTest
    NAME "Test_Tutorial_ScalingTest"
    ARGS ""
    COMM serial mpi
    NUM_MPI_PROCS 4
  )

  TRIBITS_ADD_EXECUTABLE(
    Tutorial_laplace2d
    SOURCES laplace2d.cpp
    COMM serial mpi
    )
<<<<<<< HEAD
=======
    
  set(XML_FILES_TEST ${XML_FILES_ALL})
  list(REMOVE_ITEM XML_FILES_TEST
    challenges/Elasticity2D_reference.xml
    challenges/Elasticity2D_parameters.xml
    challenges/Recirc2D_reference.xml
    challenges/Recirc2D_parameters.xml
    blocked_rebalancing.xml
    ml_ParameterList.xml
    ml_RebalancedParameterList.xml
    muelu_ParameterList.xml
    n2_easy_agg.xml
    n2_easy_agg2.xml
    n2_easy_agg3.xml
    n2_easy_export.xml
    s4av.xml
    s4c.xml
    s5a.xml
    s5b.xml
    # xmls never mentioned in the tutorial
    n1_easy_pg.xml
    s2_easy.xml
    s2b.xml
    s3a_easy.xml
    s3b_easy.xml
    s3c.xml
    s3c_easy.xml
  )
    
  foreach(xml_file IN LISTS XML_FILES_TEST)
    # without directory and without extension
    get_filename_component(xml_we "${xml_file}" NAME_WE)
>>>>>>> 77498d6d5b0 (MueLu: Tutorial: Fix challenge mode in script, using the same general TutorialDriver.cpp)

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s1_easy"
    ARGS "--xml=s1_easy.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )
<<<<<<< HEAD
=======
  endforeach()
  
  # foreach() are split so that tests are ordered by matrix type
  foreach(xml_file IN LISTS XML_FILES_TEST)
    get_filename_component(xml_we "${xml_file}" NAME_WE)
>>>>>>> 77498d6d5b0 (MueLu: Tutorial: Fix challenge mode in script, using the same general TutorialDriver.cpp)

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s1_easy_3levels_smoothed"
    ARGS "--xml=s1_easy_3levels_smoothed.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

<<<<<<< HEAD
  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s1_easy_3levels_unsmoothed"
    ARGS "--xml=s1_easy_3levels_unsmoothed.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s1_easy_exercise"
    ARGS "--xml=s1_easy_exercise.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s1_easy_jacobi"
    ARGS "--xml=s1_easy_jacobi.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s1_easy_jacobi2"
    ARGS "--xml=s1_easy_jacobi2.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s2a"
    ARGS "--xml=s2a.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s2_adv_a"
    ARGS "--xml=s2_adv_a.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s2_adv_b"
    ARGS "--xml=s2_adv_b.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s2_adv_c"
    ARGS "--xml=s2_adv_c.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s2b"
    ARGS "--xml=s2b.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s2_easy"
    ARGS "--xml=s2_easy.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3a"
    ARGS "--xml=s3a.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3a_easy"
    ARGS "--xml=s3a_easy.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3b"
    ARGS "--xml=s3b.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3b1"
    ARGS "--xml=s3b1.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3b2"
    ARGS "--xml=s3b2.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3b3"
    ARGS "--xml=s3b3.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3b_easy"
    ARGS "--xml=s3b_easy.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s3c"
    ARGS "--xml=s3c.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  # TRIBITS_ADD_TEST(
  #   Tutorial_laplace2d
  #   NAME "Test_Tutorial_laplace2d_s3c_easy"
  #   ARGS "--xml=s3c_easy.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 2
  #   )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s4a"
    ARGS "--xml=s4a.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  # TRIBITS_ADD_TEST(
  #   Tutorial_laplace2d
  #   NAME "Test_Tutorial_laplace2d_s4av"
  #   ARGS "--xml=s4av.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 2
  #   )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s4b"
    ARGS "--xml=s4b.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  # TRIBITS_ADD_TEST(
  #   Tutorial_laplace2d
  #   NAME "Test_Tutorial_laplace2d_s4c"
  #   ARGS "--xml=s4c.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 2
  #   )

  # TRIBITS_ADD_TEST(
  #   Tutorial_laplace2d
  #   NAME "Test_Tutorial_laplace2d_s5a"
  #   ARGS "--xml=s5a.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 2
  #   )

  # TRIBITS_ADD_TEST(
  #   Tutorial_laplace2d
  #   NAME "Test_Tutorial_laplace2d_s5b"
  #   ARGS "--xml=s5b.xml"
  #   COMM serial mpi
  #   NUM_MPI_PROCS 2
  #   )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s6_dep"
    ARGS "--xml=s6_dep.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

  TRIBITS_ADD_TEST(
    Tutorial_laplace2d
    NAME "Test_Tutorial_laplace2d_s6_export"
    ARGS "--xml=s6_export.xml"
    COMM serial mpi
    NUM_MPI_PROCS 2
    )

=======
>>>>>>> 77498d6d5b0 (MueLu: Tutorial: Fix challenge mode in script, using the same general TutorialDriver.cpp)
ENDIF()

# IF (${PACKAGE_NAME}_ENABLE_Epetra)

#   TRIBITS_ADD_EXECUTABLE(
#     Tutorial_recirc2d
#     SOURCES recirc2d.cpp
#     COMM serial mpi
#     )

#     TRIBITS_ADD_TEST(
#     Tutorial_recirc2d
#     NAME "Test_Tutorial_recirc2d"
#     ARGS "--xml=n1_easy.xml --xml=n1_easy_pg.xml --xml=n2_easy_agg.xml --xml=n2_easy_agg2.xml --xml=n2_easy_agg3.xml --xml=n2_easy_export.xml"
#     COMM serial mpi
#     NUM_MPI_PROCS 4
#     )

#   #TRIBITS_ADD_EXECUTABLE(
#   #  Tutorial_ScalingTestParamList
#   #  SOURCES ScalingTestParamList.cpp
#   #  COMM serial mpi
#   #  )

#   #TRIBITS_ADD_TEST(
#   #  Tutorial_ScalingTestParamList
#   #  NAME "Test_Tutorial_ScalingTestParamList"
#   #  ARGS ""
#   #  COMM serial mpi
#   #  NUM_MPI_PROCS 4
#   #  )

# ENDIF()
