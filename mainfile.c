#include<stdio.h>
#include<conio.h>
#include<math.h>

int main() {
	
  int i,j,k,degree,n,NumPair;
    
	printf("\nPlease Enter number of observation you want to enter:\n");       
    
	scanf("%d",&NumPair);
    
	double x[NumPair],y[NumPair];
    
	printf("\nEnter the values of Height:\n");                
    
	for (i=0;i<NumPair;i++){
    	scanf("%lf",&x[i]);
	}

   printf("\nEnter the values of Weight:\n");                
    
	for (i=0;i<NumPair;i++){
        scanf("%lf",&y[i]);
	}
	
	degree=2;
    
	double X[2*degree+1];                        //Array that will store the values of sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    
	for (i=0;i<2*degree+1;i++) {
        X[i]=0;
        for (j=0;j<NumPair;j++)
            X[i]=X[i]+pow(x[j],i);        //consecutive positions of the array will store N,sigma(xi),sigma(xi^2),sigma(xi^3)....sigma(xi^2n)
    }
    
    double B[degree+1][degree+2],a[degree+1];            //B is the Normal matrix(augmented) that will store the equations, 'a' is for value of the final coefficients

    for (i=0 ; i<=degree ; i++)
        for ( j=0 ; j<=degree ; j++)
            B[i][j]=X[i+j];            //Build the Normal matrix by storing the corresponding coefficients at the right positions except the last column of the matrix

    double Y[degree+1];                    
    
    for (i=0;i<degree+1;i++) {    
	    Y[i]=0;
		for (j=0;j<NumPair;j++)
            Y[i]=Y[i]+pow(x[j],i)*y[j];        //consecutive positions will store sigma(yi),sigma(xi*yi),sigma(xi^2*yi)...sigma(xi^n*yi)
	}
    
    for (i=0;i<=degree;i++)
        B[i][degree+1]=Y[i];                //load the values of Y as the last column of B(Normal Matrix but augmented)
    
	degree = degree+1;                //n is made n+1 because the Gaussian Elimination part below was for n equations, but here n is the degree of polynomial and for n degree we get n+1 equations
    
    for (i=0;i<degree;i++)                    //From now Gaussian Elimination starts(can be ignored) to solve the set of linear equations (Pivotisation)
        for (k=i+1;k<degree;k++)
            if (B[i][i]<B[k][i])
                for (j=0;j<=degree;j++) {	
                    double temp=B[i][j];
                    B[i][j]=B[k][j];
                    B[k][j]=temp;
				}
				
    for (i=0;i<degree-1;i++)            //loop to perform the gauss elimination
        for (k=i+1;k<degree;k++) {
                double t=B[k][i]/B[i][i];
                for (j=0;j<=degree;j++)
                    B[k][j]=B[k][j]-t*B[i][j];    //make the elements below the pivot elements equal to zero or elimnate the variables
            }
            
    for (i=degree-1;i>=0;i--) {                //back-substitution
	   a[i]=B[i][degree];                   //make the variable to be calculated equal to the rhs of the last equation
        for (j=0;j<degree;j++)
            if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
                a[i]=a[i]-B[i][j]*a[j];
        a[i]=a[i]/B[i][i];       
		}
   
    printf("\nThe values of the coefficients are as follows:\n");
   
     printf("cofficient of x^0 is:%lf",a[0]);
     printf("\n\n");
     printf("cofficient of x^1 is:%lf",a[1]);
     printf("\n\n");
	   printf("cofficient of x^2 is:%lf",a[2]);
     printf("\n\n");
     
	   printf("enter the number to predict: \n");
    
	   int num=0;
	
	   scanf("%d",&num);
	
     double predi = (num*num)*a[2]+num*a[1]+a[0];
    
     printf("predication value will be:%lf ",predi);
    
	return 0;
}
