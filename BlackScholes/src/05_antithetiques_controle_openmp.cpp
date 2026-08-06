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

    Result resultat;

    Statistics globalStatistics;


    #pragma omp parallel
    {

        std::random_device rd;
        std::mt19937 generator(rd());
        std::normal_distribution<double> normal(0.0,1.0);

        Statistics localStatistics;

        #pragma omp for
        for (std::size_t i = 1; i < Nsim+1; i++)
        {
            double z = normal(generator);
            
            
            double ST1 = S0*std::exp((r-0.5*sigma*sigma)*T + sigma*std::sqrt(T)*z);
            double ST2 = S0*std::exp((r-0.5*sigma*sigma)*T - sigma*std::sqrt(T)*z);

            double discountST = std::exp(-r*T)*0.5*(ST1+ST2);

            double discountPayoff = std::exp(-r*T)*0.5*(std::max(ST1-K,0.0)+std::max(ST2-K,0.0));

            updateStatistics(localStatistics,discountST,discountPayoff);
        }

        #pragma omp critical
        {
            globalStatistics = mergeStatistics(globalStatistics,localStatistics);
        }
    }


    double meanDiscountST = globalStatistics.meanX;
    double meanDiscountPayoff = globalStatistics.meanY;

    double varDiscountPayoff =  globalStatistics.M2Y/static_cast<double>(Nsim-1);
    double varDiscountST = globalStatistics.M2X/static_cast<double>(Nsim-1);
    double covDiscountPayoffDiscountST = globalStatistics.CXY/static_cast<double>(Nsim-1);

    double beta = covDiscountPayoffDiscountST/varDiscountST;

    resultat.price = meanDiscountPayoff -beta*(meanDiscountST - S0);

    double varControlDiscountPayoff = varDiscountPayoff - covDiscountPayoffDiscountST*covDiscountPayoffDiscountST/varDiscountST;

    double stdPrice = std::sqrt(varControlDiscountPayoff/static_cast<double>(Nsim));

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
        std::size_t Nsim = 1e9;

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



