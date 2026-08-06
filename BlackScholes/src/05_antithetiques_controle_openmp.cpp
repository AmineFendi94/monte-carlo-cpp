#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <limits>
#include <omp.h>
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

struct Statistics
{
    std::size_t n = 0;

    double meanX = 0.0;
    double meanY = 0.0;
    
    double M2X = 0.0;
    double M2Y = 0.0;

    double CXY = 0.0;
};


void updateStatistics(Statistics& s, double X , double Y)
{
    s.n++;
    double deltaX = X-s.meanX;
    double deltaY = Y - s.meanY;

    double n = static_cast<double>(s.n);

    s.meanX += deltaX/n;
    s.meanY += deltaY/n;

    s.M2X += deltaX*(X-s.meanX);
    s.M2Y += deltaY*(Y-s.meanY);

    s.CXY += deltaX*(Y-s.meanY);
}


Statistics mergeStatistics(Statistics const& a , Statistics const& b)
{
    if (a.n == 0) return b;
    if (b.n == 0) return a;

    Statistics s;

    s.n = a.n + b.n;

    double deltaX = b.meanX - a.meanX;
    double deltaY = b.meanY - a.meanY;

    double nB = static_cast<double>(b.n);
    double nA = static_cast<double>(a.n);
    double n = static_cast<double>(s.n);

    s.meanX = a.meanX + deltaX*nB/n;
    s.meanY = a.meanY + deltaY*nB/n;

    s.M2X = a.M2X + b.M2X + deltaX*deltaX*nA*nB/n;
    s.M2Y = a.M2Y + b.M2Y + deltaY*deltaY*nA*nB/n;

    s.CXY = a.CXY + b.CXY + deltaX*deltaY*nA*nB/n;

    return s;
}




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

    Statistics globalStatistics;

    const double drift = (r - 0.5*sigma*sigma)*T;
    const double diffusion = sigma*std::sqrt(T);
    const double discount = std::exp(-r*T);

    #pragma omp parallel
    {
        const unsigned int threadId = static_cast<unsigned int>(omp_get_thread_num());
        std::seed_seq seed{1234u, threadId};
        std::mt19937 generator(seed);
        std::normal_distribution<double> normal(0.0,1.0);

        Statistics localStatistics;

        #pragma omp for schedule(static) nowait
        for (std::size_t i = 0; i < Nsim; ++i)
        {
            const double z = normal(generator);
            const double ST1 = S0*std::exp(drift + diffusion*z);
            const double ST2 = S0*std::exp(drift - diffusion*z);
            const double discountST = discount*0.5*(ST1+ST2);
            const double discountPayoff = discount*0.5*(std::max(ST1-K,0.0)+std::max(ST2-K,0.0));

            updateStatistics(localStatistics,discountST,discountPayoff);
        }

        #pragma omp critical
        {
            globalStatistics = mergeStatistics(globalStatistics,localStatistics);
        }
    }

    if (globalStatistics.n != Nsim)
    {
        throw std::runtime_error("Le nombre de simulations fusionnées est incorrect");
    }

    const double meanDiscountST = globalStatistics.meanX;
    const double meanDiscountPayoff = globalStatistics.meanY;
    const double sampleCount = static_cast<double>(globalStatistics.n);

    const double varDiscountPayoff = globalStatistics.M2Y/(sampleCount - 1.0);
    const double varDiscountST = globalStatistics.M2X/(sampleCount - 1.0);
    const double covDiscountPayoffDiscountST = globalStatistics.CXY/(sampleCount - 1.0);

    if (varDiscountST <= std::numeric_limits<double>::epsilon())
    {
        throw std::runtime_error("La variance de la variable de contrôle est nulle");
    }

    const double beta = covDiscountPayoffDiscountST/varDiscountST;

    resultat.price = meanDiscountPayoff -beta*(meanDiscountST - S0);

    const double varControlDiscountPayoff = std::max(
        varDiscountPayoff - covDiscountPayoffDiscountST*covDiscountPayoffDiscountST/varDiscountST,
        0.0);

    const double stdPrice = std::sqrt(varControlDiscountPayoff/sampleCount);

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


