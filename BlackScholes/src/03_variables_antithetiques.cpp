#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>


struct Result
{
    double price = 0.0;
    double CILow = 0.0;
    double CIHigh = 0.0;
    double length = 0.0;
};


Result callBlackScholes(double S0, double K , double r , double T , double sigma , std::size_t Nsim)
{
    if (Nsim < 2)
    {
        throw std::invalid_argument("Nsim doit etre strictement supérieur à 1");
    }
    if (!std::isfinite(S0) || !std::isfinite(K) || !std::isfinite(r)
        || !std::isfinite(T) || !std::isfinite(sigma)
        || S0 <= 0.0 || K < 0.0 || T < 0.0 || sigma < 0.0)
    {
        throw std::invalid_argument("Les paramètres du modèle sont invalides");
    }

    Result resultat;

    std::mt19937 generator(1234u);
    std::normal_distribution<double> normal(0.0,1.0);

    double M2 = 0.0;

    const double drift = (r - 0.5*sigma*sigma)*T;
    const double diffusion = sigma*std::sqrt(T);
    const double discount = std::exp(-r*T);

    for (std::size_t i = 0; i < Nsim; ++i)
    {
        const double z = normal(generator);
        const double ST1 = S0*std::exp(drift + diffusion*z);
        const double ST2 = S0*std::exp(drift - diffusion*z);
        const double discountPayoff = discount*0.5*(std::max(ST1-K,0.0)+std::max(ST2-K,0.0));
        const double delta = discountPayoff - resultat.price;
        const double count = static_cast<double>(i + 1);

        resultat.price += delta/count;
        M2+= delta*(discountPayoff-resultat.price);
    }

    const double varDiscountPayoff = M2/static_cast<double>(Nsim-1);
    const double stdPrice = std::sqrt(varDiscountPayoff/static_cast<double>(Nsim));

    resultat.CIHigh = resultat.price + 1.96*stdPrice;
    resultat.CILow = resultat.price - 1.96*stdPrice;
    resultat.length = resultat.CIHigh - resultat.CILow;

    return resultat;
    
}

std::size_t readSimulationCount(int argc, char* argv[])
{
    constexpr std::size_t defaultCount = 1'000'000;
    if (argc == 1) return defaultCount;
    if (argc != 2) throw std::invalid_argument("Utilisation : programme [Nsim]");

    const std::string argument(argv[1]);
    if (argument.empty() || argument.front() == '-')
        throw std::invalid_argument("Nsim doit être un entier supérieur à 1");

    std::size_t parsedCharacters = 0;
    const unsigned long long value = std::stoull(argument, &parsedCharacters);
    if (parsedCharacters != argument.size() || value < 2
        || value > std::numeric_limits<std::size_t>::max())
        throw std::invalid_argument("Nsim doit être un entier supérieur à 1");
    return static_cast<std::size_t>(value);
}

int main(int argc, char* argv[])
{
    try
    {
        double S0 = 197.2;
        double r = 0.0;
        double sigma = 0.156582142;
        double K = 200;
        double T = 1.0;
        const std::size_t Nsim = readSimulationCount(argc, argv);

        auto resulat = callBlackScholes(S0, K,r,T,sigma,Nsim);


        std::cout << "Prix = " << resulat.price << '\n';
        std::cout << "CI Low = " << resulat.CILow << '\n';
        std::cout << "CI High = " << resulat.CIHigh << '\n';
        std::cout << "length = " << resulat.length << '\n';
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Erreur : " << e.what() << '\n';
        return 1;
    }
    


}


