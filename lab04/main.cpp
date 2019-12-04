#include <iostream>
#include <fstream>
#include <math.h>

#define EPS 1
#define DELTA 0.1
#define NX 150
#define NY 100
#define V1 10
#define V2 0
#define X_MAX (DELTA*NX)
#define Y_MAX (DELTA*NY)
#define SIGMA_X (0.1*X_MAX)
#define SIGMA_Y (0.1*Y_MAX)
#define TOL 1e-8


double p1(double x, double y)
    { return exp( -1 * pow( x - 0.35 * X_MAX, 2) / pow(SIGMA_X, 2) - pow( y - 0.5 * Y_MAX, 2 ) / pow(SIGMA_Y, 2) ); }

double p2(double x, double y)
    { return - exp( -1 * pow( x - 0.65 * X_MAX, 2) / pow(SIGMA_X, 2) - pow( y - 0.5 * Y_MAX, 2 ) / pow(SIGMA_Y, 2) ); }

double p(int x, int y)
    { return p1(x*DELTA,y*DELTA) + p2(x*DELTA,y*DELTA); }

void relaksacja_globalna( double omega_g ){

    double vs[NX+1][NY+1];
    double vn[NX+1][NY+1];
    double s = 0.;
    double sp = 0.;
    int counter = 0;

    std::ofstream fs;
    fs.open( "dane1.txt");


    for( int i=0; i<=NX; i++ ){
        for( int j=0; j<=NY; j++ ){
            vs[i][j] = 0.;
            vn[i][j] = 0.;
        }
    }

    for( int i=0; i<=NX; i++ ){
        vs[i][0] = 10.;
        vn[i][0] = 10.;
    }


    do{
        sp = s;

        // (9)
        for( int i=1; i<=NX-1; i++ ){
            for( int j=1; j<=NY-1; j++ )
                vn[i][j] = 1/4. * ( vs[i+1][j] + vs[i-1][j] + vs[i][j+1] + vs[i][j-1] + DELTA * DELTA / EPS * p(i,j));
        }

        // (10, 11)
        for( int j=1; j<NY; j++ ){
            vn[0][j] = vn[1][j];
            vn[NX][j] = vn[NX-1][j];
        }

        // (12)
        for( int i=0; i<=NX; i++ ){
            for( int j=0; j<=NY; j++ ){
                vs[i][j] = (1 - omega_g ) * vs[i][j] + omega_g * vn[i][j];
            }
        }

        // (17)
        s = 0;
        for( int i=0; i<=NX-1; i++ ){
            for( int j=0; j<=NY-1; j++ )
                s += pow(DELTA, 2) * ( 1/2. * pow( ( vs[i+1][j] - vs[i][j] ) / DELTA, 2 ) +  1/2. * pow( ( vs[i][j+1] - vs[i][j] ) / DELTA , 2 ) - p(i,j) * vs[i][j] );
        }
        fs << counter << " " << s << std::endl;

        counter++;
    } while( fabs( (s - sp) / sp ) > TOL );


    fs.close();
    fs.open( "dane2.txt");

    double blad[NX+1][NY+1];

    for( int i=1; i<=NX-1; i++ ){
        for( int j=1; j<=NY-1; j++ ){
            blad[i][j] = ( vn[i+1][j] - 2*vn[i][j] + vn[i-1][j] ) / pow(DELTA, 2) + ( vn[i][j+1] - 2*vn[i][j] + vn[i][j-1] ) / pow(DELTA, 2) + p(i,j)/EPS;
            std::cout << blad[i][j] << " ";
        }
        std::cout << "" << std::endl;
    }

    fs.close();

}


void relaksacja_lokalna( double omega_l ){

    double v[NX+1][NY+1];
    double s = 0.;
    double sp = 0.;
    int counter = 0;

    std::ofstream fs;
    fs.open( "dane1.txt");

    for( int i=0; i<=NX; i++ ){
        for( int j=0; j<=NY; j++ )
            v[i][j] = 0.;
    }

    for( int i=0; i<=NX; i++ )
        v[i][0] = 10.;
    
    do{
        sp = s;

        for( int i=1; i<=NX-1; i++ ){           
            for( int j=1; j<=NY-1; j++ )
                v[i][j] = (1 - omega_l) * v[i][j] + omega_l/4. * ( v[i+1][j] + v[i-1][j] + v[i][j+1] + v[i][j-1] + pow(DELTA, 2)/EPS * p(i,j));
        }

        for( int j=1; j<NY; j++ ){
            v[0][j] = v[1][j];
            v[NX][j] = v[NX-1][j];
        }

        s = 0;
        for( int i=0; i<=NX-1; i++ ){
            for( int j=0; j<=NY-1; j++ )
                s += pow(DELTA, 2) * ( 1/2. * pow( ( v[i+1][j] - v[i][j] ) / DELTA, 2 ) +  1/2. * pow( ( v[i][j+1] - v[i][j] ) / DELTA , 2 ) - p(i,j) * v[i][j] );
        }
    fs << counter << " " << s << std::endl;
        counter++;

    } while( fabs((s - sp) / sp) > TOL );


    fs.close();
    fs.open( "dane2.txt");

    double blad[NX+1][NY+1];

    for( int i=1; i<=NX-1; i++ ){
        for( int j=1; j<=NY-1; j++ ){
            blad[i][j] = ( v[i+1][j] - 2*v[i][j] + v[i-1][j] ) / pow(DELTA, 2) + ( v[i][j+1] - 2*v[i][j] + v[i][j-1] ) / pow(DELTA, 2) + p(i,j)/EPS;
            fs << blad[i][j] << " ";
        }
        fs << "" << std::endl;
    }

    fs.close();

}

int main(void){

    relaksacja_globalna( 0.6 );
    // relaksacja_globalna( 1 );

    // relaksacja_lokalna(1);
    // relaksacja_lokalna(1.4);
    // relaksacja_lokalna(1.8);
    // relaksacja_lokalna(1.9);


}
