#include <iostream>
#include <fstream>
#include <cmath>

#define NX 150
#define NT 1000
#define DELTA 0.1
#define DELTA_T 0.05
#define X_A 7.5
#define SIGMA 0.5
#define X_F 2.5


constexpr double x( int i );
constexpr double t( int j );
constexpr int delta_kroneckera(double x, double xf);
constexpr double aF(double x, double t);

constexpr void warunek_poczatkowy_u(double u[]);

void rownanie_falowe( double alpha, double beta, const char* nazwa_pliku );



int main (void){
    rownanie_falowe(0, 0  , "dane1.dat");
    rownanie_falowe(0, 0.1, "dane2.dat");
    rownanie_falowe(0, 1  , "dane3.dat");
    rownanie_falowe(1, 1  , "dane4.dat");
}



constexpr double x( int i )
    { return DELTA * i; }

constexpr double t( int j )
    { return DELTA_T*j; }

constexpr int delta_kroneckera(double x, double xf)
    { return ( fabs(x-xf) < 1e-6 ) ? 1 : 0; }

constexpr double aF(double x, double t)
    { return cos( 50.*t / (NT*DELTA_T) ) * delta_kroneckera(x, X_F); }


constexpr void warunek_poczatkowy_u(double u[]){
    for( int i=1; i<=NX-1; i++)
        u[i] = exp( -pow(x(i)-X_A, 2) / (2*pow(SIGMA, 2) ) );
}


void rownanie_falowe( double alpha, double beta, const char* nazwa_pliku ){
    double u0[NX+1]{0};
    double u [NX+1]{0};
    double v [NX+1]{0};
    double vp[NX+1]{0};
    double a [NX+1]{0};

    warunek_poczatkowy_u(u);
    for(int i=1; i<=NX-1; i++)
        u0[i] = u[i];

    for(int i=1; i<=NX-1; i++){
        a[i] = (u[i+1] - 2*u[i] + u[i-1]) / pow(DELTA, 2) - beta * ( u[i] - u0[i] ) / DELTA_T + alpha*aF(x(i), 0);
    }


    std::ofstream fs;
    fs.open(nazwa_pliku);

    for(int n=1; n<=NT; n++){
        for(int i=1; i<=NX-1; i++){
            vp[i] = v[i] + DELTA_T/2.*a[i];
            u0[i] = u[i];
            u[i] = u[i] + DELTA_T * vp[i];
        }
        for(int i=1; i<=NX-1; i++){
            a[i] = ( u[i+1] - 2*u[i] + u[i-1] ) / pow(DELTA, 2) - beta * ( u[i] - u0[i] ) / DELTA_T + alpha*aF( x(i), t(n) );        
            v[i] = vp[i] + DELTA_T/2.*a[i];
        }

        double E = DELTA/4. * ( pow( (u[1] - u[0])/DELTA , 2 ) + pow( (u[NX] - u[NX-1])/DELTA, 2 ) );
        for(int i=1; i<=NX-1; i++){
            E+= DELTA/2. * ( pow(v[i], 2) + pow( ( ( u[i+1] - u[i-1] )/( 2*DELTA ) ) , 2 ) );
        }
        // fs << t(n) << " " << E << std::endl;
        for(int i=1; i<=NX-1; i++)
            fs << t(n) << " " << x(i) << " " << u[i] << std::endl;
         fs << std::endl;
    }

}
