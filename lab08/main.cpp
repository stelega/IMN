#include <fstream>
#include <iostream>
#include <math.h>

#define DELTA 0.01
#define NX 400
#define NY 90
#define i1 200
#define i2 210
#define j1 50
#define SIGMA (10*DELTA)
#define XA 0.45
#define YA 0.45
#define IT_MAX 3500


double x(int i);
double y(int j);
void adwekcja_dyfuzja( double D );
void wypelnianie_u0( double u0[NX+1][NY+1] );
void wypelnianie_v( double vx[NX+1][NY+1], double vy[NX+1][NY+1], double psi[NX+1][NY+1] );
void wypelnianie_psi( double psi[NX+1][NY+1] );
void iteracja_Picarda( double u0[NX+1][NY+1], double u1[NX+1][NY+1], double delta_t, double D, double vx[NX+1][NY+1], double vy[NX+1][NY+1]  );


int main(void){
    adwekcja_dyfuzja(0);
    // adwekcja_dyfuzja(0.1);
}



double x(int i)
    { return i*DELTA; }


double y(int j)
    { return j*DELTA; }


void adwekcja_dyfuzja(double D){
    double u0 [NX+1][NY+1];
    double u1 [NX+1][NY+1];
    double vx [NX+1][NY+1];
    double vy [NX+1][NY+1];
    double psi[NX+1][NY+1];

    wypelnianie_psi(psi);
    wypelnianie_v(vx, vy, psi);

    // obliczanie v_max
    double v_max = 0;
    for( int i=0; i<=NX; i++ ){
        for( int j=0; j<=NY; j++ ){
            if( sqrt( pow( vx[i][j], 2 ) + pow( vy[i][j], 2 ) ) > v_max )
                v_max = sqrt( pow( vx[i][j], 2 ) + pow( vy[i][j], 2 ) );
        }
    }

    double delta_t = DELTA/(4.*v_max);
    wypelnianie_u0(u0);

    iteracja_Picarda(u0, u1, delta_t, D, vx, vy);

}


void wypelnianie_u0( double u0[NX+1][NY+1] ){
    for( int i=0; i<=NX; i++ ){
        for( int j=0; j<=NY; j++ )
            u0[i][j] = 1./( 2*M_PI*pow(SIGMA, 2) ) * exp( - ( pow( x(i)-XA, 2 ) + pow( y(j)-YA, 2 ) ) / ( 2*pow(SIGMA, 2) ) ) ;
    }
}


void wypelnianie_v( double vx[NX+1][NY+1], double vy[NX+1][NY+1], double psi[NX+1][NY+1] ){
    for( int i=1; i<=NX-1; i++ ){
        for( int j=1; j<=NY-1; j++ ){
            vx[i][j] =  ( psi[i][j+1] - psi[i][j-1] ) / ( 2.*DELTA );
            vy[i][j] = -( psi[i+1][j] - psi[i-1][j] ) / ( 2.*DELTA );
        }
    }
    for( int i=i1; i<=i2; i++ ){
        for( int j=0; j<=j1; j++ ){
            vx[i][j] = 0.;
            vy[i][j] = 0.;
        }
    }
    for( int i=1; i<=NX-1; i++ ){
        vx[i][0 ] = 0;
        vx[i][NY] = 0;
        vy[i][0 ] = 0;
        vy[i][NY] = 0;       
    }
    for( int j=0; j<=NY; j++ ){
        vx[0][j] = vx[1][j];
        vx[NX][j] = vx[NX-1][j];
    }

}


void wypelnianie_psi( double psi[NX+1][NY+1] ){

    std::ifstream fs1;
    fs1.open("psi.dat");

    if(!fs1)
        std::cout << "Nie znaleziono pliku" << std::endl;
    else{
        int i, j;

        while (fs1 >> i){
            fs1 >> j;
            fs1 >> psi[i][j];
        }
    }

    fs1.close();
}


void iteracja_Picarda( double u0[NX+1][NY+1], double u1[NX+1][NY+1], double delta_t, double D,  double vx[NX+1][NY+1], double vy[NX+1][NY+1] ){

    std::ofstream fs2, fs3;
    fs2.open("dane.dat");


    for( int it=1; it<=IT_MAX; it++ ){
        for( int i=0; i<=NX; i++ ){
            for( int j=0; j<=NY; j++ ){
                u1[i][j]=u0[i][j];
            }
        }

        for( int k=1; k<=20; k++ ){
            for( int i=0 ;i<=NX; i++ ){
                for( int j=1 ;j<=NY-1; j++){
                    if( i<i1 || i>i2 || j>j1 ){
                        if(i==0){
                            u1[i][j] = ( 1./( 1+( 2*D*delta_t / pow(DELTA, 2)) ) ) * ( u0[i][j] - (delta_t/2.) * vx[i][j] *
                            ( ( (u0[i+1][j] - u0[NX][j])/(2.*DELTA) ) + (u1[i+1][j] - u1[NX][j])/(2.*DELTA) ) - (delta_t / 2) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.*DELTA) + (u1[i][j+1] - u1[i][j-1])/(2.*DELTA) ) + delta_t/2. * D * 
                            ( ( u0[i+1][j] + u0[NX][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(DELTA,2) + ( u1[i+1][j] + u1[NX][j] + u1[i][j+1] + u1[i][j-1] )/pow(DELTA,2) )
                            );
                        }
                        else if(i==NX){
                            u1[i][j] = ( 1./( 1+( 2*D*delta_t / pow(DELTA, 2)) ) ) * ( u0[i][j] - (delta_t/2.) * vx[i][j] *
                            ( ( (u0[0][j] - u0[i-1][j])/(2.*DELTA) ) + (u1[0][j] - u1[i-1][j])/(2.*DELTA) ) - (delta_t / 2) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.*DELTA) + (u1[i][j+1] - u1[i][j-1])/(2.*DELTA) ) + delta_t/2. * D * 
                            ( ( u0[0][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(DELTA,2) + ( u1[0][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(DELTA,2) )
                            );
                        }
                        else{
                            u1[i][j] = ( 1./( 1+( 2*D*delta_t / pow(DELTA, 2)) ) ) * ( u0[i][j] - (delta_t/2.) * vx[i][j] *
                            ( ( (u0[i+1][j] - u0[i-1][j])/(2.*DELTA) ) + (u1[i+1][j] - u1[i-1][j])/(2.*DELTA) ) - (delta_t / 2) * vy[i][j] * 
                            ( ( u0[i][j+1] - u0[i][j-1] )/(2.*DELTA) + (u1[i][j+1] - u1[i][j-1])/(2.*DELTA) ) + delta_t/2. * D * 
                            ( ( u0[i+1][j] + u0[i-1][j] + u0[i][j+1] + u0[i][j-1] - 4*u0[i][j] )/pow(DELTA,2) + ( u1[i+1][j] + u1[i-1][j] + u1[i][j+1] + u1[i][j-1] )/pow(DELTA,2) )
                            );
                        }
                    }
                }
            }
        }        
        
        for( int i=0; i<=NX; i++ ){
            for( int j=0; j<=NY; j++ )
                u0[i][j]=u1[i][j];
        }
        double c=0.;
        double xsr=0.;
        for( int i=0; i<=NX; i++ ){
            for( int j=0; j<=NY; j++ ){
                c += u0[i][j];
                xsr += x(i) * u0[i][j];
            }
        }
        fs2 << it*delta_t << " " << c* pow(DELTA, 2) << " " << xsr*pow(DELTA, 2) << std::endl;

    }

    fs2.close();
}
