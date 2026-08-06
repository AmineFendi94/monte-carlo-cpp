#include <iostream>
#include <cmath>
#include <random>


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


    for (std::size_t i = 1; i < Nsim+1; i++)
    {
        double z = normal(generator);
        
        double ST = S0*std::exp((r-0.5*sigma*sigma)*T + sigma*std::sqrt(T)*z);

        double discountPayoff = std::exp(-r*T)*std::max(ST-K,0.0);

        resultat.price += discountPayoff;
    }

    resultat.price /= static_cast<double>(Nsim);

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


        std::cout << "Price = " << resulat.price << std::endl;
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error : " << e.what() << '\n';
        return 1;
    }
    


}



