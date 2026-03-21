#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>


using namespace std;

//Problem Set 2:
// Excersice 1 | Question 1 | i)

//Pricing a Europian Put Option with the Binomial Method
//Input variables of the corresponding pricing function are : 
// S,K,r,T,σ, number of time steps N

double EuropeanPut_Binomial(double S, double K, double r,double v,double T,int N){
	// δt=T/N;
	double deltaT = T/N ;
	// u = exp(σ*(sqrt(δt))
	double u = exp( v* sqrt(deltaT) );
	// d = 1/u 
	double d = 1/u ;
	// p = (e^(rδt)-d)/(u-d) 
	double p = ( exp(r*(deltaT))-d ) / (u-d) ;
	
	/*Here is the key aspect of the Lattice. We have the paths of up and down movements 
	 for the pricing of the European Put so in order to include any single one we need 
	 to store these prices as a data structure of an Array !*/
	 
	 double lattice[N+1][N+1]; 
	 
	 // now we need a for loop to obtain the payoffs based on i states and N steps :
	 // This means we obtain the payoff at the maturity == Expiration Date == Excersice time
	 
	 // f_{iN} = max(0,S_0 * u^i*d^(N-1) - K) where this function represents every possible payment at the end
	 
	 for(int i = 0; i<N+1; i++){
	 	lattice[i][N] = max(K- S* pow(u,i)* pow(d,N-i),0.0 );
	 }
	 //We constructed the above for loop because the procedure of obtaining the prices is done backwards in time
	 //Starting from the last node = maturity and by the Present Value we are going to compute the possible payoff
	 //At the time and step -> (time,step) == (i,j) nodes we desire. 
	 
	 //For that we need another for loop :
	 for( int j=N-1; j>-1; j--){
	 	for(int i=0; i<j+1; i++){
	 		lattice[i][j] = exp(-r*deltaT)* ( p* lattice[i+1][j+1] + (1-p)*lattice[i][j+1]);
		 }
	 }
	 return lattice[0][0]; //This is the price at the t=0 -> starting time 
	
}

int main(){
	int N = 5;
	double S = 35;
	double K = 40;
	double r = 0.03;
	double v = 0.25;
	double T = 0.75;
	
	double put = EuropeanPut_Binomial(S,K,r,v,T,N);
	
	cout << "The price of the European Put via the Binomial Method is: "<< put ;
	
	return 0;
}

