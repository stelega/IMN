#include <fstream>
#include <iostream>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#define DELTA 1
#define NX 40
#define NY 40
#define N ((NX+1)*(NY+1))
#define DELTA_T 1
#define TA 40
#define TB  0
#define TC 30
#define TD  0
#define KB 0.1
#define KD 0.6
#define IT_MAX 2000
#define l ( i + j*(NX+1) )


double x(int i);
double y(int j);
void wypelnianie_a( gsl_matrix *a);
void wypelnianie_b( gsl_matrix *b);
void wypelnianie_c( gsl_vector *c);
void wypelnianie_T( gsl_vector *T);
double nabla_kwadrat_razy_T( gsl_vector *T, int k );


int main(void){

    gsl_matrix *a = gsl_matrix_calloc(N, N);
    gsl_matrix *b = gsl_matrix_calloc(N, N);
    gsl_vector *c = gsl_vector_calloc(N);
    gsl_vector *d = gsl_vector_calloc(N);
    gsl_vector *T = gsl_vector_calloc(N);
    gsl_permutation* p = gsl_permutation_calloc(N);
    int signum = 0;

    wypelnianie_a(a);
    wypelnianie_b(b);
    wypelnianie_c(c);
    wypelnianie_T(T);


    gsl_linalg_LU_decomp(a, p, &signum);

    std::ofstream f1, f2;
    // tu program wpisuje x, y oraz T(x, y)
    f1.open("dane1.dat");
    // tu program wpisuje x, y oraz pow(nabla, 2)*T(x, y)
    f2.open("dane2.dat");

    for( int it=1;  it<=IT_MAX; it++ ){
        gsl_blas_dgemv( CblasNoTrans, 1, b, T, 0, d );
        gsl_blas_daxpy( 1, c, d );
        gsl_linalg_LU_solve( a, p, d, T );
        if( it == 100 || it == 200 || it == 500 || it == 1000 || it == 2000 ){
            std::cout << it << std::endl;
            for( int i=1; i<=NX-1; i++ ){
                for( int j=1; j<=NY-1; j++ ){
                    f1 << x(i) << " " << y(j) << " " << gsl_vector_get(T, l) << std::endl;
                    f2 << x(i) << " " << y(j) << " " << nabla_kwadrat_razy_T(T, l) << std::endl;
                }
                f1 << std::endl;
                f2 << std::endl;
            }
            f1 << std::endl;
            f1 << std::endl;
            f2 << std::endl;
            f2 << std::endl;
        }
    }

    f1.close();
    f2.close();

    gsl_matrix_free(a);
    gsl_matrix_free(b);
    gsl_vector_free(c);
    gsl_vector_free(T);
    gsl_permutation_free(p);
}



double x(int i)
    { return i*DELTA; }


double y(int j)
    { return j*DELTA; }


void wypelnianie_a( gsl_matrix *a){
    // srodek
    for( int i=1; i<=NX-1; i++){
        for(int j=1; j<=NY-1; j++){
            gsl_matrix_set(a, l, l-NX-1, DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(a, l, l-1,    DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(a, l, l+1,    DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(a, l, l+NX+1, DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(a, l, l,      -(2*DELTA_T)/pow(DELTA,2) - 1);
        }
    }

    // lewy i prawy brzeg
    int i = 0;
    for(int j=0; j<=NY; j++)
        gsl_matrix_set(a, l, l, 1);
    i = NX;
    for(int j=0; j<=NY; j++)
        gsl_matrix_set(a, l, l, 1);

    // gorny brzeg
    int j = NY;
    for( i=1; i<=NX-1; i++ ){
        gsl_matrix_set(a, l, l-NX-1, -1 / (KB * DELTA));
        gsl_matrix_set(a, l, l, 1 + 1 / (KB * DELTA));
    }

    // dolny brzeg
    j = 0;
    for( i=1; i<=NX-1; i++ ){
        gsl_matrix_set(a, l, l,      1 + 1 / (KD * DELTA));
        gsl_matrix_set(a, l, l+NX+1, -1 / (KD * DELTA));
    }

}


void wypelnianie_b( gsl_matrix *b){
    // srodek
    for( int i=1; i<=NX-1; i++){
        for(int j=1; j<=NY-1; j++){
            gsl_matrix_set(b, l, l-NX-1, -DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(b, l, l-1,    -DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(b, l, l+1,    -DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(b, l, l+NX+1, -DELTA_T / (2*pow(DELTA, 2)));
            gsl_matrix_set(b, l, l,      (2*DELTA_T)/pow(DELTA,2) - 1);
        }
    }

    // lewy i prawy brzeg
    int i = 0;
    for(int j=0; j<=NY; j++)
        gsl_matrix_set(b, l, l, 1);

    i = NX;
    for(int j=0; j<=NY; j++)
        gsl_matrix_set(b, l, l, 1);

    // gorny brzeg
    int j = NY;
    for( int k=0; k<N; k++ ){
        for( i=1; i<=NX-1; i++ )
        gsl_matrix_set(b, l, k, 0);
    }
    
    // dolny brzeg
    j = 0;
    for( int k=0; k<N; k++ ){
        for( i=1; i<=NX-1; i++ )
        gsl_matrix_set(b, l, k, 0);
    }

}


void wypelnianie_c( gsl_vector *c){
    int i=0;
    for(int j=0; j<=NY; j++)
        gsl_vector_set(c, l, 0);
        
    i=NX;
    for(int j=0; j<=NY; j++)
        gsl_vector_set(c, l, 0);
    
    int j=NY;
    for(int i=1; i<=NX-1; i++)
        gsl_vector_set(c, l, TB);

    j=0;
    for(int i=1; i<=NX-1; i++)
        gsl_vector_set(c, l, TD);

}

void wypelnianie_T( gsl_vector *T){
    int i=0;
    for(int j=0; j<=NY; j++)
        gsl_vector_set(T, l, TA);
   
    i=NX;
    for(int j=0; j<=NY; j++)
        gsl_vector_set(T, l, TC);

    for(i=1; i<=NX-1; i++){
        for(int j=0; j<=NY; j++)
        gsl_vector_set(T, l, 0);
    }
    
}


double nabla_kwadrat_razy_T( gsl_vector *T, int k ){
    return ( ( gsl_vector_get(T, k+1) - 2 * gsl_vector_get(T, k) + gsl_vector_get(T, k-1)  )/pow(DELTA, 2))
    + ( ( gsl_vector_get(T, k+NX+1) - 2*gsl_vector_get(T, k) + gsl_vector_get(T, k-NX-1) )/pow(DELTA, 2) );
}

