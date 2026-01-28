#include <iostream>
#include "roots.hpp"
#include <cmath>

bool bisection(std::function<double(double)> f,
               double a, double b,
               double *root){
                const double tol = 1e-7;
                const double max_iterations = 1000000;
                double fa= f(a);
                double fb= f(b);

                if (fa * fb > 0) {   
                    return false;
                }
                
                double c=0;
                int iterations = 0;

                while(iterations < max_iterations){
                    c=(a+b)/2.0;

                    double fc = f(c);

                    if ( fabs(fc) < tol ){
                        *root = c;
                        return true;
                    }
                    else if (fa * fc < 0){
                        b=c;
                        fb=fc;
                    } 
                    else {
                        a = c;
                        fa = fc;
                    } 
                    iterations++;  
                }
                std::cout <<"Max iteration reached ("<<max_iterations<<"))\n";
                *root = (a+b)/2.0;
                return false;   
               }



bool regula_falsi(std::function<double(double)> f,
                  double a, double b,
                  double *root){

                const double tol = 1e-7;
                const double max_iterations = 1000000;
                double fa= f(a);
                double fb= f(b);

                 if (fa * fb > 0) {  //want opposite signs 
                    return false;
                }

                double c;
                int iterations = 0;

                while(iterations < max_iterations){
                    c = a-(fa*(b-a))/(fb-fa);
                    double fc = f(c);

                    if ( fabs(fc) < tol ){
                        *root = c;
                        return true;
                    }
                    else if (fa * fc < 0){
                        b=c;
                        fb=fc;
                    } 
                    else {
                        a = c;
                        fa = fc;
                    } 
                    iterations++;  
                }


                return false;
                  }

bool newton_raphson(std::function<double(double)> f,
                    std::function<double(double)> g,
                    double a, double b, double c,
                    double *root) {
                        const double tol = 1e-6;
                        const double max_iterations = 1000;

                        double x=c;
                        double x_new;

                        for ( int i =0; i< max_iterations; i++){
                            double gx=g(x);
                            if (fabs(gx) < 1e-12) {
                                return false; 
                            }

                             x_new = x - f(x)/gx;
                            
                            if(fabs(x_new-x)<=tol){
                                *root = x_new;
                                return true;
                            }
                            x = x_new;
                            
                        }
                        
                        return false;
                    
                }


bool secant(std::function<double(double)> f,
            double a, double b, double c,
            double *root) {
                const int max_iterations = 1000;
                const double tol = 1e-6;

                if (c < a || c > b) {
                    return false; 
                }
                double x0 = c;
                double x1 = c + 0.1;

                if (x1 < a || x1 > b){
                    return false; 
                    }
                
                for (int i = 0; i < max_iterations; i++) {
                    if(f(x1)-f(x0) == 0){
                        return false; 
                    }

                    double x2 = x1 - f(x1) * (x1 - x0) / (f(x1) - f(x0));

                    if (fabs(x2 - x1) < tol) {
                        *root = x2;
                        return true;
                    }

                    x0=x1;
                    x1=x2;
                }
            

                return false;
            }



               



