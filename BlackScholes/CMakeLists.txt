add_executable(black_scholes src/BlackScholes.cpp)
add_executable(black_scholes_interval_confiance src/BlackScholesIntervalConfiance.cpp)
add_executable(black_scholes_antithetic src/BlackScholesAntithetic.cpp)
add_executable(black_scholes_control_variate src/BlackScholesVC.cpp)
add_executable(black_scholes_openmp src/BlackScholesOpenMP.cpp)

foreach(target
        black_scholes
        black_scholes_interval_confiance
        black_scholes_antithetic
        black_scholes_control_variate
        black_scholes_openmp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
endforeach()

find_package(OpenMP REQUIRED)
target_link_libraries(black_scholes_openmp PRIVATE OpenMP::OpenMP_CXX)
