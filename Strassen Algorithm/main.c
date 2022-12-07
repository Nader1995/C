
// *******************************************************************************
// ***                                                                         *** 
// ***                                                                         ***
// ***       Strassen algorithm, Version # 0.30                                ***
// ***       Rectangular matrix multiplication mod 2                           ***
// ***                                                                         ***
// ***       Here I develope strassen algorithm                                ***
// ***       Nader Safari                                                      ***
// ***       Date: 13991221, 20210311                                          ***
// ***       Copy Left 2021, N. Safari, All Left Reserved!                     ***
// ***                                                                         *** 
// *******************************************************************************


#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// ===============================================================================
// functions

// ================================================
// power ( finds p in 2^p > n)
						
int power(int base, int exp){
    
	int i=0, p=1;
    
	while(++i<exp+1 && 1*(p*=base));
    return p;
}

// ================================================
// matrix adder 

void Madd(int* A, int* B, int* C, int n, int x){
	
    int i,j, m= x>0?n/2:n;
    
    for(i= 0; i<m; i++){
    	
        for(j= 0; j<m; j++)
            *(C+i*m+j) = *(A+i*n+j) + *(B+i*n+j);
	}
}

// ================================================
// matrix subtracter 

void Msub(int* A, int* B, int* C, int n, int x){
	
    int i,j, m= x>0?n/2:n;
    
	for(i= 0; i<m; i++){
     
	    for(j= 0; j<m; j++){
         
		    *(C+i*m+j) = *(A+i*n+j) - *(B+i*n+j);
		}
	}
}

// ================================================
// random generator	

int randGenerator(int lower, int upper){
	
	int num = 0;	
	num = (rand() % (upper - lower + 1)) + lower; 
	return num;
}

// ================================================
// strassen algorithm

void strassen(int* A, int* B, int* C, int n){
	
    int i, j, k;
    
	if(n==2){	
		
		int P = 0, Q = 0, R = 0, S = 0, T = 0, U = 0, V = 0;
		
		if(*A != 0){
			
			R=(*(B+1)-*(B+n+1));
		}
		
		if(*(A+n+1) != 0){
			
			S=(*(B+n)-*B);
		}
		
		if(*B != 0){
			
			Q=(*(A+n)+*(A+n+1));
		}
		
		if(*(B+n+1) != 0){
			
			T=(*A+*(A+1));
		}
		
        P=(*A+*(A+n+1))*(*B+*(B+n+1));  //P=(A[0][0]+A[1][1])*(B[0][0]+B[1][1])
//        Q=(*(A+n)+*(A+n+1))*(*B);   //Q=(A[1][0]+A[1][1])*B[0][0]
//        R=(*A)*(*(B+1)-*(B+n+1));   //R=A[0][0]*(B[0][1]-B[1][1])
//        S=(*(A+n+1))*(*(B+n)-*B);   //S=A[1][1]*(B[1][0]-B[0][0])
//        T=(*A+*(A+1))*(*(B+n+1));   //T=(A[0][0]+A[0][1])*B[1][1]
        U=(*(A+n)-*A)*(*B+*(B+1));  //U=(A[1][0]-A[0][0])*(B[0][0]+B[0][1])
        V=(*(A+1)-*(A+n+1))*(*(B+n)+*(B+n+1));  //V=(A[0][1]-A[1][1])*(B[1][0]+B[1][1])
		
		
		
        *C=(P+S-T+V);
        *(C+1)=(R+T);
        *(C+n)=(Q+S);
        *(C+n+1)=(P+R-Q+U);
        
        /*
        int P1=0, Q1=0, P2=0, Q2=0, P3=0, Q3=0, P4=0, Q4=0;
        
        if(*A != 0){
        	 P1 = *B;
        	 Q1 = *(B+1);
		}
		
		if(*(A+1) != 0){
        	 P2 = *(B+n);
        	 Q2 = *(B+n+1);
		}
		
		if(*(A+n) != 0){
        	 P3 = *B;
        	 Q3 = *(B+1);
		}
		printf("\n%d\n", *(A+n+1));
		if(*(A+n+1) != 0){
        	 P4 = *(B+n);
        	 Q4 = *(B+n+1);
		}
		
        *C= P1 + P2;
        *(C+1)= Q1 + Q2;
        *(C+n)= P3 + P4;
        *(C+n+1)= Q3 + Q4;
        */
    }

    else{

        int m=n/2, x[m][m], y[m][m], o[n][n];

        for(i= 0; i<n; i++){
        	
        	for(j= 0; j<n; j++){
            	
            	o[i][j]= 0;
			}
		}
                
        /*P=(A[0][0]+A[1][1])*(B[0][0]+B[1][1])*/
        
		int P[m][m];
        
		Madd(A, A+m*(n+1), x, n, 1);
        Madd(B, B+m*(n+1), y, n, 1);
        strassen(x, y, P, m);

        /*Q=(A[1][0]+A[1][1])*B[0][0]*/
        
		int Q[m][m];
        
		Madd(A+m*n, A+m*(n+1), x, n, 1);
        Madd(B, o, y, n, 1);
        strassen(x, y, Q, m);

        /*R=A[0][0]*(B[0][1]-B[1][1])*/
        
		int R[m][m];
        
		Madd(A, o, x, n, 1);
        Msub(B+m, B+m*(n+1), y, n, 1);
        strassen(x, y, R, m);

        /*S=A[1][1]*(B[1][0]-B[0][0])*/
        
		int S[m][m];
        
		Madd(A+m*(n+1), o, x, n, 1);
        Msub(B+m*n, B, y, n, 1);
        strassen(x, y, S, m);

        /*T=(A[0][0]+A[0][1])*B[1][1]*/
        
		int T[m][m];
        
		Madd(A, A+m, x, n, 1);
        Madd(B+m*(n+1), o, y, n, 1);
        strassen(x, y, T, m);

        /*U=(A[1][0]-A[0][0])*(B[0][0]+B[0][1])*/
        
		int U[m][m];
        
		Msub(A+m*n, A, x, n, 1);
        Madd(B, B+m, y, n, 1);
        strassen(x, y, U, m);

        /*V=(A[0][1]-A[1][1])*(B[1][0]+B[1][1])*/
        
		int V[m][m];
        
		Msub(A+m, A+m*(n+1), x, n, 1);
        Madd(B+m*n, B+m*(n+1), y, n, 1);
        strassen(x, y, V, m);


        /*Calculating the 4 parts for the result matrix*/
        
		int W[m][m], X[m][m], Y[m][m], Z[m][m];
        
		Msub(V,T,x,m,0);
        Madd(S,x,y,m,0);
        Madd(P,y,W,m,0); // W=P+S-T+V
        Madd(R,T,X,m,0); // X==R+T
        Madd(Q,S,Y,m,0); // Y=Q+S
        Msub(U,Q,x,m,0);
        Madd(R,x,y,m,0);
        Madd(P,y,Z,m,0); // Z=P+R-Q+U

        /*Conquering 4 parts in the result matrix*/
        
		for (i=0;i<m;i++){
			
            for (j=0;j<m;j++){
                *(C+i*n+j) = W[i][j]; //C[0][0]=W
                *(C+i*n+j+m) = X[i][j]; //C[0][1]=X
                *(C+(i+m)*n+j) = Y[i][j]; //C[1][0]=Y
                *(C+(i+m)*n+j+m) = Z[i][j]; //C[1][1]=Z
            }
		}
    }
}

// ===============================================================================
// main part

void main(){
	
	// ==============================================================
	// Include random seed
	
	srand(time(0));
	
	// ==============================================================
	// Variables & pointers
	
    int i, j, k, m, n, n1=300 , n2=320 , n3=400 , n4, o= 0;
	int lower, upper;
    int A[n1][n2];
    int B[n2][n3];
    int mul[n1][n3];
    
    // ===============================================================
	// Here I have to produce random number 0 and 1 ...
	// ... and put them in matrices a and b
	
	lower= 0;
	upper= 1;
	for(i= 0; i<n1; i++){
		
		for(j= 0; j<n2; j++){
			
			A[i][j] = randGenerator(lower, upper);
		}
	}
	    
	for(i= 0; i<n2; i++){
		
		for(j= 0; j<n3; j++){
			
			B[i][j] = randGenerator(lower, upper);
		}
	}
	
	// if condition n1>n2 is true return n1, otherwise return n2
  	
    n4=n1>n2?n1:n2;
    
    // likewise
  	
	n=n3>n4?n3:n4;
	
	while(n>(m=power(2,++o)));
	
	int a[m][m], b[m][m], C[m][m];
    
    // build the square 2^m * 2^m matrices
    
    for(i= 0; i<m; i++){
    	
		for(j= 0; j<m; j++){
            a[i][j]= 0;
            b[i][j]= 0;
        }
	}
		
    for(i= 0; i<n1; i++){
    	
    	for(j= 0; j<n2; j++){
    		
    	a[i][j]= A[i][j];	
		}
   	}
        
    for(i= 0; i<n2; i++){
    	
    	for(j= 0; j<n3; j++){
        	
        	b[i][j]= B[i][j];
		}
	}
	
	
    // ===============================================================
	// call the main function
	
	strassen(a,b,C,m);    //Calling the function.
	
	// ===============================================================
	// tester
	
	for(i=0; i<n1; i++){
		
		for(j=0; j<n3; j++){
			
			mul[i][j] = 0;
			
			for(k=0; k<n2; k++){
				
				mul[i][j] += A[i][k] * B[k][j];
			}
			mul[i][j] = mul[i][j] % 2; 
		}
	}

    
	/*Printing test matrix (Decimal method)*/
    printf("\nThis is the first matrix:");
    for(i=0;i<n1;i++){
        printf("\n\n\n");
        for(j=0;j<n3;j++)
            printf("\t%d",mul[i][j]);
    }

	/*Printing the final matrix*/
    printf("\n\n\nThis is the final matrix:");
    for(i=0;i<n1;i++){
        printf("\n\n\n");
        for(j=0;j<n3;j++)
            printf("\t%d",C[i][j]);
    }    
        
}


