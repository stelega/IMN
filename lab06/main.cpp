#include <iostream>
#include <fstream>
#include <math.h>
#include "mgmres.c"
#include "mgmres.h"

#define DELTA 0.1
#define TOL 1e-8

typedef double ( *P_PTR )(double,double,double,double);


double p1( double x, double y, double x_max, double y_max){
    double sigma = x_max/10;
    return exp( pow((x-0.25*x_max) ,2)/ pow(sigma,2) - pow((y-0.5*y_max) ,2)/ pow(sigma,2)  );    
}

double p2( double x, double y, double x_max, double y_max ){
    double sigma = x_max/10;
    return exp( pow((x-0.75*x_max) ,2)/ pow(sigma,2) - pow((y-0.5*y_max) ,2)/ pow(sigma,2)  );    
}

double p0( double x, double y, double x_max, double y_max){
    return 0;
}


void fun( int nx, int ny, double eps1, double eps2, double v1, double v2, double v3, double v4, P_PTR p1, P_PTR p2 ){



    int N = (nx+1)*(ny+1);
    double V[N];
    double a[5*N];
    double b[N];
    int ia[N+1];
    int ja[5*N];

    double y_max = ny*DELTA;
    double x_max = nx*DELTA;  
    double eps[N];
    for (int l = 0; l <= N + 1; l++){
        int j = floor(l / (nx + 1));
        int i = l - j * (nx + 1);
        if (i <= nx / 2)
            eps[l] = eps1;
        else
            eps[l] = eps2;
    }
    
    for(int i=0; i<N+1; i++)
        ia[i] = -1;

    for(int i=0; i<5*N; i++){
        a[i] = 0;
        ja[i] = 0;
    }

    for(int i=0; i<N; i++){
        V[i] = 0;
        b[i] = 0;
    }

    int k = -1;
    for(int l=0; l<N ;l++){
        int j = floor( l/(nx+1) );
        int i = l-j*(nx+1);
        int brzeg=0;   // wskaźnik  położenia: 0-środek  obszaru; 1-brzeg
        double vb=0.; // potencjal  na  brzegu
        if(i==0){ //lewy  brzeg
            brzeg=1;
            vb=v1;
        }
        if(j==ny){ //górny  brzeg
        brzeg=1;
        vb=v2;
        }
        if(i==nx){ //prawy  brzeg
        brzeg=1;
        vb=v3;
        }
        if(j==0) { //dolny  brzeg
        brzeg=1;
        vb=v4;
        }
        // wypełniamy  od razu  wektor  wyrazów  wolnych
        b[l]= - (p1(i*DELTA, j*DELTA, x_max, y_max) + p2(i*DELTA, j*DELTA, x_max, y_max)); //jeśli w środku  jest  gęstość
        if(brzeg ==1) {
            b[l]=vb; // wymuszamy  wartość  potencjału  na  brzegu
        }
    
        ia[l]=-1; // wskaźnik  dla  pierwszego  el. w wierszu
    //lewa  skrajna  przekatna
        if(l-nx-1>=0 &&  brzeg ==0 ){
            k++;
            if(ia[l]<0)
                ia[l]=k;
            a[k]=eps[l+1] / pow(DELTA, 2) ;;
            ja[k]=l-nx-1;
        }
        // poddiagonala
        if(l-1>=0 &&  brzeg ==0 ){
            k++;
            if(ia[l]<0)
                ia[l]=k;
            a[k]=eps[l] / pow(DELTA, 2) ;;
            ja[k]=l-1;
        }
        // diagonala
            k++;
        if(ia[l]<0)
            ia[l]=k;
        if(brzeg ==0)
            a[k]= -( 2*eps[l] + eps[l+1] + eps[l+nx+1] ) / pow(DELTA, 2) ;
        else
            a[k]=1;

        ja[k]=l;
        // naddiagonala
        if(l<N &&  brzeg ==0 ){
            k++;
            a[k]=eps[l+1] / pow(DELTA, 2) ;
            ja[k]=l+1;
        }
        // prawa  skrajna  przekątna
        if(l<N-nx-1 && brzeg ==0 ){
            k++;
            a[k]= eps[l+nx+1] / pow(DELTA, 2) ;
            ja[k]=l+nx+ 1;
        }
    }

    int nz_num=k+1; //ilosc  niezerowych  elementow  (1  element  ma  indeks  0)
    ia[N]= nz_num;
        
    double itr_max= 500;
    double mr = 500;
    double tol_abs = 1e-8;
    double tol_rel= 1e-8;
    pmgmres_ilu_cr( N, nz_num, ia, ja, a, V, b, itr_max, mr, tol_abs, tol_rel  );


    std::ofstream fs1, fs2, fs3;
    fs1.open("dane1.txt");


    for(int i=0; i<=nx; i++){
        for(int j=0; j<=ny; j++)
            fs1 << i*(nx+1)+j << " " << i << " " << j << " " << b[i*(nx+1)+j] << std::endl;
    }


    fs2.open("dane2.txt");

    for(int i=0; i<N*5; i++){
        fs2 << i << " " << a[i] << std::endl;
    }
    
    fs3.open("dane3.txt");

    for(int i=0; i<nx; i++){
        for(int j=0; j<ny; j++){
            int l = i + j*(nx+1);
            fs3 << i << " " << j << " " << V[l] << std::endl;
        }
    }

    fs1.close();
    fs2.close();
    fs3.close();
    

}



int main(void){

    // fun( 4, 4, 1, 1, 10, -10, 10, -10, p0, p0 );
    fun( 50, 50, 1, 1, 10, -10, 10, -10, p0, p0 );
    // fun( 200, 200, 1, 1, 10, -10, 10, -10, p0, p0 );
    // fun( 100, 100, 1, 1, 10, -10, 10, -10, p0, p0 );

    // fun( 100, 100, 1, 1, 0, 0, 0, 0, p1, p2 );
    
    

}
