#include <iostream>
#include <stdio.h>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <random>

using namespace std;

//Problem Set 2 : Question 3 | i)

/*We are tasked to compute an exotic Call Option named Barrier Knockout 
with the technique of Monte Carlo. First we will address the key 
aspects of this type of Option.

1) For the Barrier Knockout type we know that its payoff depends on wether
or not the underlying's price will reach or exceed a certain predetermined price
In the case of Knock-Out => It is active for excersice the whole time 
up to the moment where the price reaches or exceeds that barrier. Once there it "terminates"

2) Here because the "life" of the Option depends on a threshold we will need it's value.
So we name it Sb == Barrier value in order to compare it with S to see if we can or not 
excersice it. !

We will again use the uniform distribution alongside the inverse of normal cdf.
*/

//We need this number to be Uniform 

double unifRand(){
	return rand()/double(RAND_MAX); 
	/*Here we scale him properly to be always in [0,1] space*/
}

/* We take now the Inverse CDF of the Normal distribution N(μ,σ^2) 
that we need in order to generate the Normal Random Numbers that we need 
for the MC simulation */
double inverse_of_normal_cdf(const double p, const double mu, const double sigma)
{
    if (p <= 0.0 || p >= 1.0)
    {
        std::stringstream os;
        os << "Invalid input argument (" << p
            << "); must be larger than 0 but less than 1.";
        throw std::invalid_argument(os.str());
    }

    double r, val;

    const double q = p - 0.5;

    if (std::abs(q) <= .425) {
        r = .180625 - q * q;
        val =
            q * (((((((r * 2509.0809287301226727 +
                33430.575583588128105) * r + 67265.770927008700853) * r +
                45921.953931549871457) * r + 13731.693765509461125) * r +
                1971.5909503065514427) * r + 133.14166789178437745) * r +
                3.387132872796366608)
            / (((((((r * 5226.495278852854561 +
                28729.085735721942674) * r + 39307.89580009271061) * r +
                21213.794301586595867) * r + 5394.1960214247511077) * r +
                687.1870074920579083) * r + 42.313330701600911252) * r + 1);
    }
    else {
        if (q > 0) {
            r = 1 - p;
        }
        else {
            r = p;
        }

        r = std::sqrt(-std::log(r));

        if (r <= 5) 
        {
            r += -1.6;
            val = (((((((r * 7.7454501427834140764e-4 +
                .0227238449892691845833) * r + .24178072517745061177) *
                r + 1.27045825245236838258) * r +
                3.64784832476320460504) * r + 5.7694972214606914055) *
                r + 4.6303378461565452959) * r +
                1.42343711074968357734)
                / (((((((r *
                    1.05075007164441684324e-9 + 5.475938084995344946e-4) *
                    r + .0151986665636164571966) * r +
                    .14810397642748007459) * r + .68976733498510000455) *
                    r + 1.6763848301838038494) * r +
                    2.05319162663775882187) * r + 1);
        }
        else { /* very close to  0 or 1 */
            r += -5;
            val = (((((((r * 2.01033439929228813265e-7 +
                2.71155556874348757815e-5) * r +
                .0012426609473880784386) * r + .026532189526576123093) *
                r + .29656057182850489123) * r +
                1.7848265399172913358) * r + 5.4637849111641143699) *
                r + 6.6579046435011037772)
                / (((((((r *
                    2.04426310338993978564e-15 + 1.4215117583164458887e-7) *
                    r + 1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r + .0148753612908506148525)
                    * r + .13692988092273580531) * r +
                    .59983220655588793769) * r + 1);
        }

        if (q < 0.0) {
            val = -val;
        }
    }

    return mu + sigma * val;
}

//Now we generate X_i random var. from the N(0,1)

double Stand_Normal_Rand(){
	double u = unifRand(); //Step 1 : Gen u_i from Unif(0,1)
	
	while(u<=0 || u>=1){
		u= unifRand();
	}
	return inverse_of_normal_cdf(u,0,1);
}

//We proceed by creating the Function of the Barrier Knockout Put Option 

//We take as base the code provided for the Barier KnockIn Call Price.


double sample_var = 0.0;

//Pricing the Barrier Knockout Put Option via the Monte Carlo method :

double BarrierKnockoutPut_MC(double S,double Sb,double K,double r,double v,double T,int N,int num_sims){
	
	double dt = T/N ;
	double nudt = ( r - 0.5*v*v )*dt;
	double vdt = v*sqrt(dt);
	double disc_payoff_sum = 0.0 ;
	double disc_payoff_squared_sum = 0.0;
	
	double SPath[N+1];
	
	int i;
	for(i=0; i<num_sims; i++){
		SPath[0] = S;
		int Knocked = 0;
		double disc_payoff = 0.0; 
		int j;
		for(j=0; j<N; j++){
			double epsilon = Stand_Normal_Rand();
			SPath[j+1] = SPath[j] * exp(nudt + vdt*epsilon);
			if( Sb > S && SPath[j+1] >= Sb ){ // This is the case where we have an Up and Out Barrier 
		       Knocked = 1;
		}else if(Sb < S && SPath[j+1] <= Sb){ //This is the case of the Down and Out Barrier
	           Knocked = 1;
	}
}
	if(Knocked==0){
		disc_payoff = exp(-r*T)*max(K- SPath[N],0.0); // Put Option Payoff 
	}
	double disc_payoff_squared = disc_payoff * disc_payoff; 
	disc_payoff_sum += disc_payoff;
	disc_payoff_squared_sum += disc_payoff_squared;
		
}
    double disc_payoff_average = disc_payoff_sum / num_sims; 
    sample_var = (disc_payoff_squared_sum - num_sims*disc_payoff_average*disc_payoff_average)/(num_sims -1);
    return disc_payoff_average; 

}


int main(){
	//Setting up the values of the actual parameters 
	
	double S= 40;
	double Sb = 55;
	double K = 39; 
	double r = 0.05;
	double T = 0.5;
	double v = 0.3;
	int N = 1000;
	int num_sims = 150000;
	
	double price = BarrierKnockoutPut_MC(S,Sb,K,r,v,T,N,num_sims);
	
	cout << "The price of the Barrier Knockout Put via the Monte Carlo is: " << price << endl;
	
	return 0 ;
	
}