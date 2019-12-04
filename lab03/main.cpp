#include <fstream>
#include <iostream>
#include "math.h"
#include <algorithm>

#define T0 0
#define X0 0.01
#define V0 0
#define DELTA_T0 1
#define S 0.75
#define P 2
#define T_MAX 40
#define ALPHA 5
#define M_DELTA 1e-10



double f(double t, double x, double v){
    return v;
}


double g(double t, double x, double v){
    return ALPHA*( 1 - x*x) * v - x;
}

double F( double xn, double xn1, double vn, double vn1, double delta_t ){
    return xn1 - xn - delta_t/2. * ( f( 1., xn, vn ) + f( 1., xn1, vn1 ) );
}

double G( double xn, double xn1, double vn, double vn1, double delta_t ){
    return vn1 - vn - delta_t/2. * ( g( 1., xn, vn ) + g( 1., xn1, vn1 ) );
}

void metoda_trapezow( double xn, double vn, double delta_t, double* xn1, double* vn1 ){
    
    double delta_x, delta_v;    
    double a11 = 1;
    double a12 = (-1)*delta_t/2.;

    double xk = xn;
    double vk = vn;

    do{
        double a21 = (-1)*delta_t/2. * ( -2 * ALPHA  * xk * vk - 1 );
        double a22 = 1 - delta_t/2. * ALPHA *( 1 - xk*xk  );

        delta_x = ( (-1)* F(xn, xk, vn, vk, delta_t) * a22 - (-1) * G(xn, xk, vn, vk, delta_t) * a12 )/ ( a11*a22 - a12*a21 );
        delta_v = ( a11* (-1) * G(xn, xk, vn, vk, delta_t) - a21 * (-1) * F(xn, xk, vn, vk, delta_t) )/ ( a11*a22 - a12*a21 );
        xk += delta_x;
        vk += delta_v;
    }
    while( fabs(delta_x) > 1e-10 || fabs(delta_v) > 1e-10 );

    *xn1 = xk;
    *vn1 = vk;

}


void metoda_RK2( double xn, double vn, double delta_t, double* xn1, double* vn1 ){


    double k1x, k2x, k1v, k2v;

    k1x = f( 1 , xn, vn );
    k1v = g( 1 , xn, vn );
    
    k2x = f( 1 + delta_t , xn + delta_t * k1x , vn + delta_t * k1v );
    k2v = g( 1 + delta_t , xn + delta_t * k1x , vn + delta_t * k1v );

    *xn1 = xn + delta_t/2. *( k1x + k2x );
    *vn1 = vn + delta_t/2. *( k1v + k2v );
   
}


void kontrola_kroku_czasowego_trapezy( double TOL ){
 
    double xn = X0;
    double vn = V0;
    double xn1, vn1, xn2, vn2, xn3, vn3 ;
    double delta_t = DELTA_T0;
    double Ex, Ev;
    double t = 0;

    std::ofstream fs;
    fs.open( "dane2.txt" );
    
    do{
        metoda_trapezow( xn, vn, delta_t, &xn1, &vn1 );
        metoda_trapezow( xn1, vn1, delta_t, &xn2, &vn2 );

        metoda_trapezow( xn, vn, 2*delta_t, &xn3, &vn3 );

        Ex = ( xn2 - xn3 ) / (pow(2, P) - 1.);
        Ev = ( vn2 - vn3 ) / (pow(2, P) - 1.); 

        if( std::max( fabs( Ex ), fabs(Ev) ) < TOL ){
            t += 2*delta_t;
            xn = xn2;
            vn = vn2;
            fs << t << " " << delta_t << " " << xn << " " << vn << "\n";
        }

        delta_t = pow( ( S * TOL )/ std::max( fabs( Ex ), fabs(Ev) ), 1./(P+1) ) * delta_t; 
    }
    while ( t<T_MAX );


    fs.close();

}


void kontrola_kroku_czasowego_RK2( double TOL ){

    double xn = X0;
    double vn = V0;
    double xn1, vn1, xn2, vn2, xn3, vn3 ;
    double delta_t = DELTA_T0;
    double Ex, Ev;
    double t = 0;

    std::ofstream fs;
    fs.open( "dane1.txt" );
    
    do{
        metoda_RK2( xn, vn, delta_t, &xn1, &vn1 );
        metoda_RK2( xn1, vn1, delta_t, &xn2, &vn2 );

        metoda_RK2( xn, vn, 2*delta_t, &xn3, &vn3 );

        Ex = ( xn2 - xn3 ) / (pow(2, P) - 1.);
        Ev = ( vn2 - vn3 ) / (pow(2, P) - 1.); 

        if( std::max( fabs( Ex ), fabs(Ev) ) < TOL ){
            t += 2*delta_t;
            xn = xn2;
            vn = vn2;
            fs << t << " " << delta_t << " " << xn << " " << vn << "\n";
        }

        delta_t = pow( ( S * TOL )/ std::max( fabs( Ex ), fabs(Ev) ), 1./(P+1) ) * delta_t; 
    }
    while ( t<T_MAX );


    fs.close();

}



int main(void){
    
    kontrola_kroku_czasowego_RK2(1e-2);
    // kontrola_kroku_czasowego_RK2(1e-5);
    kontrola_kroku_czasowego_trapezy(1e-2);
    // kontrola_kroku_czasowego_trapezy(1e-5);
}
