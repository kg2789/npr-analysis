set (TEST_DATA_DIR ${PROJECT_SOURCE_DIR}/../test_data/)

add_executable(test_fierz ${CMAKE_CURRENT_SOURCE_DIR}/test_fierz.cpp ${ANALYSIS_SOURCES})
target_link_libraries(test_fierz ${GRID_LIBS})
add_test(fierz_mixed_check test_fierz ${TEST_DATA_DIR}/mixed_test)

add_executable(test_freeprop ${CMAKE_CURRENT_SOURCE_DIR}/test_freeprop.cpp ${ANALYSIS_SOURCES})
target_link_libraries(test_freeprop ${GRID_LIBS})
add_test(freeprop_check test_freeprop ${TEST_DATA_DIR}/freeprop_test)

add_executable(test_jackknife ${CMAKE_CURRENT_SOURCE_DIR}/test_jackknife.cpp)
target_link_libraries(test_jackknife ${GRID_LIBS})
add_test(jackknife test_jackknife)

