#include <iostream>
#include <stdio.h>
#include <cmath>
#include <cstdlib>
#include <algorithm>

using namespace std;

//Problem Set 2:
// Excersice 1 | Question 1 | ii)

//Here we need to perform the Binomial Method for the American Put
//The main difference with the European one is that Americans Puts
//Can be excersiced prior to the maturity at any moment => t<=T. 

//So the American will have aswell the European's same inpute variables :

double AmericanPut_Binomial(double S, double K, double r,double v,double T,int N){
	// δt=T/N;
	double deltaT = T/N ;
	// u = exp(σ*(sqrt(δt))
	double u = exp( v* sqrt(deltaT) );
	// d = 1/u 
	double d = 1/u ;
	// p = (e^(rδt)-d)/(u-d) 
	double p = ( exp(r*deltaT)-d) / (u-d) ;
	 
	double lattice[N+1][N+1]; 
	
	// Again the terminal payoffs for the American One are given through the below function: 
	// f_{iN} = max(0,S_0 * u^i*d^(N-1) - K)
	 
	for(int i = 0; i<N+1; i++){
	 	lattice[i][N] = max(K- S* pow(u,i)* pow(d,N-i), 0.0 );
	}
	// The only thing we need to change is in the following for loop of the Backward Procedure:
	for(int j=N-1; j>-1; j--){
		for(int i=0; i<j+1; i++){
		/*The main 2 aspects of the American Put that we implement in this loop are:
		It's price if it wasnt excercised up until maturity and the Immediate Excercise Value
		
		With these 2 we have the American Put Option's value : max(noExcercise,ImmediateExcercise)
		*/
		
		double noExcerciseValue = exp(-r*deltaT)* ( p* lattice[i+1][j+1] + (1-p)*lattice[i][j+1]);
			
		//Underlying price at the node(i,j):
		double S_ij = S * pow(u,i) * pow(d,j-i);
		//Payoff given at t<=T <=> Immediate Excercise 
		double ImmediateExcersiceValue = max(K-S_ij, 0.0);
		
		//American Put Option based on the nodes (i,j):
		lattice[i][j] = max(noExcerciseValue,ImmediateExcersiceValue);
		}
	}
	return lattice[0][0];
}


int main(){
	
	int N = 5;
	double S = 35;
	double K = 40;
	double r = 0.03;
	double v = 0.25;
	double T = 0.75;
	
	double put = AmericanPut_Binomial(S,K,r,v,T,N);
	
	cout << "The price of the American Put via the Binomial Method is: "<< put <<endl;
	
	
	//The price of the European Call ATM:
	//To be At The Money we need the K==S 
	
	double eu_atm = EUCall(S,S,r,T,v);
	
	cout << "The price of the European Call Option in ATM state is: " << eu_atm ;
	return 0;
}
	 
	 
	 
	 