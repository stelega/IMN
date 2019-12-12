#include <fstream>
#include <iostream>
#include <math.h>

#define DELTA 0.01
#define ro 1.
#define mi 1.
#define NX 200
#define NY 90
#define i1 50
#define j1 55
#define IT_MAX 20000
#define DELTA_X DELTA
#define DELTA_Y DELTA


void relaksacja( int Qwe );
double x( int x );
double y( int y );
void wypelnianie_V( double V[NX+1][NY+1], double Qwe, double Qwy );
void wypelnianie_C( double C[NX+1][NY+1], double V[NX+1][NY+1], double Qwe, double Qwy );



int main(void){
    relaksacja(-1000);
    // relaksacja(-4000);
    // relaksacja(4000);
}



void relaksacja( int Qwe ){
    double V[NX+1][NY+1] {0};
    double C[NX+1][NY+1] {0};

    double Qwy = Qwe * ( pow( y(NY) ,3) - pow(y(j1) ,3) - 3*y(j1) * pow(y(NY), 2) + 3 * pow(y(j1), 2) * y(NY) ) / pow(y(NY), 3);

    wypelnianie_V(V, Qwe, Qwy);

    // brzeg A dla C
    for( int j=j1; j<=NY; j++)
        C[0][j] = Qwe/(2*mi) * ( 2*y(j) - y(j1) - y(NY));

    // brzeg C dla C
    for( int j=0; j<=NY; j++)
        C[NX][j] = Qwy/(2*mi) * ( 2*y(j) - y(NY) );

    std::ofstream fs1;
    fs1.open("dane.dat");
    
    int omega;
    for (int IT=1; IT<=IT_MAX; IT++){
        if(IT < 2000) 
            omega = 0;
        else
            omega = 1;
    
        wypelnianie_C(C, V, Qwe, Qwy);
        for (int i=1; i<=NX-1; i++){
            for (int j=1; j<=NY-1; j++){
                if( i>i1 || j>j1 ){
                    V[i][j]= 1./4. * ( V[i+1][j] + V[i-1][j] + V[i][j+1] + V[i][j-1] - pow(DELTA, 2) * C[i][j] ) ;
                    C[i][j]= 1./4. * ( C[i+1][j] + C[i-1][j] + C[i][j+1] + C[i][j-1] ) - omega * ro/(16*mi) *
                    ( ( V[i][j+1]-V[i][j-1] )*(C[i+1][j]-C[i-1][j])  - (V[i+1][j]-V[i-1][j])*(C[i][j+1]-C[i][j-1]) ) ;
                }
            }
        }
        double T = 0;
        for(int i=1; i<=NX-1; i++)
            T += fabs(V[i+1][j1+2] + V[i-1][j1+2] + V[i][j1+3] + V[i][j1+1] - 4*V[i][j1+2] - pow(DELTA, 2) * C[i][j1+2]);
    }
    double u[NX+1][NY+1] {0};
    double v[NX+1][NY+1] {0};

    for(int i=1; i<=NX-1; i++){
        for(int j=1; j<=NY-1; j++){
            if( i>i1 || j>j1 ){
                u[i][j] =  ( V[i][j+1] - V[i][j-1] ) / (2*DELTA);
                v[i][j] = -( V[i+1][j] - V[i-1][j] ) / (2*DELTA);
            }
        }
    }
    for(int i=0; i<=NX; i++){
        for(int j=0; j<=NY; j++)
            fs1 << x(i) << " " << y(j) << " " << u[i][j] << " " << v[i][j] << " " << V[i][j] << " " << C[i][j] << std::endl;

        fs1 << std::endl;
    }

    fs1.close();
}


double x( int x )
    { return x * DELTA_X; }


double y( int y )
    { return y * DELTA_Y; }


void wypelnianie_V( double V[NX+1][NY+1], double Qwe, double Qwy ){
    
    // brzeg A 
    for( int j=j1; j<=NY; j++)
        V[0][j] = Qwe/(2*mi) * ( pow(y(j),3)/3. - pow(y(j),2)/2. * ( y(j1) + y(NY) ) + y(j)*y(j1)*y(NY) );
    
    // brzeg C
    for( int j=0; j<=NY; j++)
        V[NX][j] = Qwy/(2*mi) * ( pow(y(j),3)/3. - pow(y(j),2)/2. * y(NY) ) + ( Qwe * pow(y(j1), 2) * ( 3*y(NY) - y(j1) ) )/(12*mi);

    // brzeg B
    for(int i=1; i<=NX-1; i++)
        V[i][NY] = V[0][NY];

    // brzeg D
    for(int i=i1; i<=NX-1; i++)
        V[i][0] = V[0][j1];

    // brzeg E
    for(int j=1; j<=j1; j++)
        V[i1][j] = V[0][j1];

    // brzeg F
    for(int i=1; i<=i1; i++)
        V[i][j1] = V[0][j1];
}


void wypelnianie_C( double C[NX+1][NY+1], double V[NX+1][NY+1], double Qwe, double Qwy ){

    // brzeg B
    for(int i=1; i<=NX-1; i++)
        C[i][NY] = 2./pow(DELTA,2) * ( V[i][NY-1] - V[i][NY] );

    // brzeg D
    for(int i=i1+1; i<=NX-1; i++)
        C[i][0] = 2./pow(DELTA,2) * ( V[i][1] - V[i][0] );

    // brzeg E
    for(int j=1; j<=j1-1; j++)
        C[i1][j] = 2./pow(DELTA,2) * ( V[i1+1][j] - V[i1][j] );

    // brzeg F
    for(int i=1; i<=i1; i++)
        C[i][j1] = 2./pow(DELTA,2) * ( V[i][j1+1] - V[i][j1] );

    // wierzcholek E/F
    C[i1][j1] = 1./2. * ( C[i1-1][j1] + C[i1][j1-1] );
}
