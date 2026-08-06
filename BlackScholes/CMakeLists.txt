add_executable(01_monte_carlo src/01_monte_carlo.cpp)
add_executable(02_intervalle_confiance_welford src/02_intervalle_confiance_welford.cpp)
add_executable(03_variables_antithetiques src/03_variables_antithetiques.cpp)
add_executable(04_variable_controle src/04_variable_controle.cpp)
add_executable(05_antithetiques_controle_openmp src/05_antithetiques_controle_openmp.cpp)

foreach(target
        01_monte_carlo
        02_intervalle_confiance_welford
        03_variables_antithetiques
        04_variable_controle
        05_antithetiques_controle_openmp)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
endforeach()

find_package(OpenMP REQUIRED)
target_link_libraries(05_antithetiques_controle_openmp PRIVATE OpenMP::OpenMP_CXX)
