#include <iostream>
#include <fstream>
#include <math.h>

#define DELTA 0.2
#define NX 128
#define NY 128
#define X_MAX (DELTA*NX)
#define Y_MAX (DELTA*NY)
#define TOL 1e-8



void relaksacja(){

    double v[NX+1][NY+1];

    for( int i=0; i<=NY; i++ ){
        v[0][i]  = sin( M_PI * i*DELTA / Y_MAX );
        v[NX][i] = sin( M_PI * i*DELTA / Y_MAX );
    }

    for( int i=0; i<=NX; i++ ){
        v[i][0]  =  sin( 2*M_PI *i*DELTA / X_MAX );
        v[i][NY] = -sin( 2*M_PI *i*DELTA / X_MAX );
    }

    for( int i=1; i<NX; i++ ){
        for( int j=1; j<NY; j++ ){
            v[i][j] = 0.;
        }
    }

    std::ofstream fs1, fs2;
    fs1.open("dane1.txt");
    fs2.open("dane2.txt");

    double s  = 0;
    double sp = 0;
    int k;

    for( k=16; k>=1; k/=2 ){
        int counter = 0;
        do{
            sp = s;
            for( int i=k; i<NX; i+=k ){
                for( int j=k; j<NY; j+=k ){
                    v[i][j] = 1./4. * ( v[i+k][j] + v[i-k][j] + v[i][j+k] + v[i][j-k] );
                }
            }
            s = 0;
            for( int i=0; i<NX; i+=k ){
                for( int j=0; j<NY; j+=k ){
                    s+= pow(k*DELTA, 2)/2. * ( pow( ( v[i+k][j] - v[i][j] ) / ( 2*k*DELTA ) + ( v[i+k][j+k] - v[i][j+k] ) / ( 2*k*DELTA ) , 2 ) +
                    pow( ( v[i][j+k] - v[i][j] ) / ( 2*k*DELTA ) + ( v[i+k][j+k] - v[i+k][j] ) / ( 2*k*DELTA ) , 2 ) );
                }
            }
            
            counter++;

            fs1 << k << " " << counter << " " << s << std::endl;
        } while( fabs( (s - sp) / sp ) > TOL  );


        if( k>1 ){
            for( int i=0; i<NX; i+=k ){
                for( int j=0; j<NY; j+=k ){
                    v[i+k/2][j+k/2] = 1./4. * ( v[i][j]   + v[i+k][j] + v[i][j+k] + v[i+k][j+k] );
                    if( i+k < NX )
                        v[i+k]  [j+k/2] = 1./2. * ( v[i+k][j] + v[i+k][j+k] );
                    if( j+k < NY )
                        v[i+k/2][j+k]   = 1./2. * ( v[i][j+k] + v[i+k][j+k] );
                }
            }
        }

        for( int i=0; i<=NX; i+=k ){
            for( int j=0; j<=NY; j+=k ){
                fs2 << i* DELTA << " " << j*DELTA << " " << v[i][j] << std::endl;
            }
            fs2 << std::endl;
        }
        fs2 << std::endl;
        


    }

    fs1.close();
    fs2.close();
    

}



int main(void){

    relaksacja();


}
