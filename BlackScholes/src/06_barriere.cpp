#include <algorithm>
#include <iostream>
#include <cmath>
#include <random>
#include <stdexcept>
#include <eigen3/Eigen/Dense>


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
    // Un thread peut ne recevoir aucune itération si Nsim est plus petit que
    // le nombre de threads. Dans ce cas, il n'y a aucune statistique à fusionner.
    if (a.n == 0)
    {
        return b;
    }
    if (b.n == 0)
    {
        return a;
    }

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


Result DownOutCallBlackScholes(double S0, double K , double L ,double r , double T , double sigma , std::size_t Nsim, std::size_t M )
{
    if (Nsim < 2)
    {
        throw std::invalid_argument("Nsim doit etre strictement supérieur à 1");
    }

    if (M < 1)
    {
        throw std::invalid_argument("M doit etre strictement positif");
    }

    Result resultat;

 
    double dt = T/M;
    double drift = (r-0.5*sigma*sigma)*dt;
    double diffusion = sigma*std::sqrt(dt);
    double discount = std::exp(-r*T);

    double logL = std::log(L);

    Statistics globalStatistics;

    


    #pragma omp parallel
    {

        std::random_device rd;
        std::mt19937 generator(rd());
        std::normal_distribution<double> normal(0.0,1.0);

        Statistics localStatistics;

        #pragma omp for
        for (size_t i = 0; i < Nsim; i++)
        {
            double logS1 = std::log(S0);
            double logS2 = std::log(S0);

            bool barriere1 = S0>L;
            bool barriere2 = S0>L;

            for (std::size_t j= 0; j < M; j++)
            {   
                double z = normal(generator);

                logS1 += drift + diffusion*z;
                logS2 += drift - diffusion*z;

                if (barriere1 && logS1 <=logL)
                {
                    barriere1 = false;
                }

                if (barriere2 && logS2 <=logL)
                {
                   barriere2 = false; 
                }

                if (!barriere1 && ! barriere2)
                {
                    double t = (j+1)*dt;

                    double zFinal = normal(generator);

                    logS1 +=  (r-0.5*sigma*sigma)*(T-t) + sigma*std::sqrt(T-t)*zFinal;

                    logS2 +=  (r-0.5*sigma*sigma)*(T-t) - sigma*std::sqrt(T-t)*zFinal;

                    break;
                }


            }

            double ST1 = std::exp(logS1);
            double ST2 = std::exp(logS2);

            double discountST = discount*0.5*(ST1+ST2);


            double discountPayoff1 = discount*(std::max(ST1-K,0.0)*barriere1);
            double discountPayoff2 = discount*(std::max(ST2-K,0.0)*barriere2);

            double discountPayoff = 0.5*(discountPayoff1+discountPayoff2);

            updateStatistics(localStatistics,discountST, discountPayoff);

        }

         #pragma omp critical
        {
            globalStatistics = mergeStatistics(globalStatistics,localStatistics);
        }
        
    }

    double n = static_cast<double>(globalStatistics.n);

    double meanDiscountST = globalStatistics.meanX;

    std::cout << "Martingal Check : " << " E[exp(-r*T)*ST/S0] = " << meanDiscountST/S0 << std::endl;

    double varDiscountST = globalStatistics.M2X/(n-1.0);
    double meanDiscountPayoff = globalStatistics.meanY;

    double covDiscountSTDiscountPayoff = globalStatistics.CXY/(n-1.0);
    double beta = covDiscountSTDiscountPayoff/varDiscountST;

    resultat.price = meanDiscountPayoff - beta*(meanDiscountST - S0);
    
    double varDiscountPayoff = globalStatistics.M2Y/(n-1.0);

    double varControlDiscountPayoff = std::max(
        varDiscountPayoff - covDiscountSTDiscountPayoff*covDiscountSTDiscountPayoff/varDiscountST,
        0.0);


    double stdPrice = std::sqrt(varControlDiscountPayoff/n);

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
        double L = 180;
        double T = 1.0;
        std::size_t Nsim = 1e6;
        std::size_t M = 252;

        auto resulat = DownOutCallBlackScholes(S0,K,L,r,T,sigma,Nsim,M);


        std::cout << "Prix = " << resulat.price << std::endl;
        std::cout << "CI Low = " << resulat.CILow << std::endl;
        std::cout << "CI High = " << resulat.CIHigh << std::endl;
        std::cout << "length = " << resulat.length << std::endl;
        return 0;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error : " << e.what() << '\n';
        return 1;
    }
    


}
