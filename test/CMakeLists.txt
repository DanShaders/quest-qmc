set(SOURCES
  config.cpp
)

include(FetchContent)

FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/1d17ea141d2c11b8917d2c7d029f1c4e2b9769b2.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(unit_tests ${SOURCES})

target_link_libraries(unit_tests PRIVATE GTest::gtest_main dqmc)

include(GoogleTest)
gtest_discover_tests(unit_tests)
