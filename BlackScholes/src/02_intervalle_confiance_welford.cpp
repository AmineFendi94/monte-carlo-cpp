#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>


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

    Result resultat;

    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<double> normal(0.0,1.0);

    double varDiscountPayoff = 0.0;

    // Ces constantes ne dépendent pas du tirage aléatoire.
    double drift = (r-0.5*sigma*sigma)*T;
    double diffusion = sigma*std::sqrt(T);
    double discount = std::exp(-r*T);

    for (std::size_t i = 1; i <= Nsim; i++)
    {
        double z = normal(generator);
        
        double ST = S0*std::exp(drift + diffusion*z);

        double discountPayoff = discount*std::max(ST-K,0.0);

        double delta = discountPayoff - resultat.price;



        resultat.price += delta/static_cast<double>(i);
        varDiscountPayoff+= delta*(discountPayoff-resultat.price);

    }

    varDiscountPayoff /= static_cast<double>(Nsim-1);

    double stdPrice = std::sqrt(varDiscountPayoff/static_cast<double>(Nsim));

    resultat.CIHigh = resultat.price + 1.96*stdPrice;
    resultat.CILow = resultat.price - 1.96*stdPrice;
    resultat.length = resultat.CIHigh - resultat.CILow;

    return resultat;
    
}


int main()
{
    try
    {
        double S0 = 197.2;
        double r = 0.0;
        double sigma = 0.156582142;
        double K = 200;
        double T = 1.0;
        std::size_t Nsim = 1e6;

        auto resulat = callBlackScholes(S0, K,r,T,sigma,Nsim);


        std::cout << "Prix = " << resulat.price << std::endl;
        std::cout << "CI Low = " << resulat.CILow << std::endl;
        std::cout << "CI High = " << resulat.CIHigh << std::endl;
        std::cout << "length = " << resulat.length << std::endl;
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Erreur : " << e.what() << '\n';
        return 1;
    }
    


}

