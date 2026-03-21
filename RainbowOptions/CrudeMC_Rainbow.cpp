#include <iostream>
#include <cmath>
#include <stdio.h>
#include <sstream>
#include <algorithm>
// Include <cstdlib> that gives the rand,srand
#include <cstdlib>

using namespace std;

/* Pricing a Rainbow Option via Monte Carlo
The Rainbow is of type : "Call on min"
The holder has the right to purchase the maximum asset at the strike price at expiry.
Call on max -> max(min(S1,S2,...S_n) - K, 0.0))

-------- Crude Monte Carlo Method ---------
In this example the Rainbow is constructed by only 2 assets -> S1,S2. 
We do not correlate them for start. This means -> ρ_1,2 = 0.
Thus the Cov(S1,S2) = 0 -> Might be independence or perhaps a non-linear relationship. 
But since these 2 assets are following the Geometric Brownian Motion then their log-returns
are jointly Gaussian. So now indeed :
Cov(S1,S2) = 0 && log-returns jointly Gaussian -> Independence */

// First lets set the uniform distribution number generator:
double UnifRand(){
	return rand() / double(RAND_MAX);
	//Scaling has been perfrom by dividing with RAND_MAX
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
	double u = UnifRand(); //Step 1 : Gen u_i from Unif(0,1)
	
	while(u<=0 || u>=1){
		u= UnifRand();
	}
	return inverse_of_normal_cdf(u,0,1);
}

double sample_var = 0.0;
// Constructing the Monte Carlo Pricing Function 
// for Rainbow 

double MC_Call_Price_Rainbow(double S1,double S2,double K,double r,double v1,double v2,double T,int num_sims){
	
	// The correlation is assumed to be ρ = 0
	// Setting up the variables of the Si(T) = Si(0)* exp( (r-0.5* v_i^2)T + v_i*sqrt(T)*ε_i)
	
	// Each asset here must have its own drift and volatility under the risk-neutral measure
	double nu_1 = ( r - 0.5 * v1 * v1) * T;
	double nu_2 = ( r - 0.5 * v2 * v2) * T;
	
	double v_1 = v1*sqrt(T);
	double v_2 = v2*sqrt(T);
    
    double disc_payoff_sum = 0.0;
    double disc_payoff_squared_sum = 0.0;
    
    // Main loop for the path prices generation
    int i;
    for( i=0; i<num_sims; i++){
    	
    	// Generating each ε_i ~ N(0,1) for asset i 
    	// Here there are only 2 needed.
    	
    	double epsilon1 = Stand_Normal_Rand();
    	double epsilon2 = Stand_Normal_Rand();
    	
    	double S1_T = S1 * exp(nu_1 + v_1 * epsilon1);
    	double S2_T = S2 * exp(nu_2 + v_2 * epsilon2);
    	
    	//Setting up the payoff of the 2 assets Rainbow Price
    	double minimum = min(S1_T,S2_T);
    	double disc_payoff = max(minimum - K,0.0) * exp( - r * T);
    	
    	disc_payoff_sum += disc_payoff ;
    	disc_payoff_squared_sum += disc_payoff * disc_payoff;   	
	}
	
	double disc_payoff_average = disc_payoff_sum / num_sims ;
	
	// Unbiased sample variance 
	sample_var = (disc_payoff_squared_sum - num_sims * disc_payoff_average * disc_payoff_average) / (num_sims -1);
	
	return disc_payoff_average ;
} 

int main() {
	/* Important Notice :
    In order to achieve reproducibility across the
	differnet methods : Crude MC, Antithetic Variates, Control Variates
	and hence reduce extreme deviations cause by the random noise introduced
	by the Random Number Generators we insert the following :
	This doesnt affect the paths : Same seed != Same Paths 
	*/
    //unsigned int seed = 123456u; // pick any fixed integer seed
	//srand(seed);
	
	// Getting the actual parameters :
	
    int num_sims = 4000000;   // Number of simulated asset paths
    double S1 = 100.0;      // Stock1 price
    double S2 = 100.0;      // Stock2 price
    double K = 100.0;      // Strike price
    double r = 0.05;       // Risk-free rate (5%)
    double v1 = 0.2;        // Volatility of the underlying (20%)
    double v2 = 0.2;        // Volatility of the underlying (20%)
    double T = 1.0;        // One year until expiration

    
    // Via the function we compute the Call price of Rainbow Option
    
    double call_price = MC_Call_Price_Rainbow(S1,S2,K,r,v1,v2,T,num_sims);
    
    // Relative error for alpha = 5% (95% confidence)
    double z = inverse_of_normal_cdf(0.975, 0.0, 1.0);
    double rel_error = z * std::sqrt(sample_var / num_sims) / call_price;
    
    // Print out the results 
    cout << "Number of Paths for Crude MC: " << num_sims << endl;
    cout << "Underlying 1: " << S1 << "\n";
    cout << "Underlying 2: " << S2 << "\n";
    cout << "Strike: " << K << "\n";
    cout << "Risk-Free Rate: " << r << "\n";
    cout << "Volatility1: " << v1 << "\n";
    cout << "Volatility2: " << v2 << "\n";
    cout << "Maturity: " << T << "\n";
    cout << "Correlation p between S1 AND S2: 0"<< "\n";

    cout << "Call Price: " << call_price << "\n";
    cout << "Relative Error: " << rel_error << "\n";

    return 0;
}