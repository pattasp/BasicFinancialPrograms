#include <iostream>
#include <cmath>
#include <stdio.h>
#include <sstream>
#include <algorithm>
#include <cstdlib>
#include <iomanip>


using namespace std;

/* Pricing a Rainbow Option via Monte Carlo
The Rainbow is of type : "Call on min"
The holder has the right to purchase the maximum asset at the strike price at expiry.
Call on max -> max(min(S1,S2,...S_n) - K, 0.0))

-------- Control Variates Monte Carlo --------------------------------------------
*/
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
// for Rainbow Option 

double MC_Call_Price_Rainbow_CV(double S1,double S2,double K,double r,double v1,double v2,double T,int num_sims,int num_pilot){
	
	// The correlation is assumed to be ρ = 0
	// Setting up the variables of the Si(T) = Si(0)* exp( (r-0.5* v_i^2)T + v_i*sqrt(T)*ε_i)
	
	// Each asset here must have its own drift and volatility under the risk-neutral measure
	double nu_1 = ( r - 0.5 * v1 * v1) * T;
	double nu_2 = ( r - 0.5 * v2 * v2) * T;
	
	double v_1 = v1*sqrt(T);
	double v_2 = v2*sqrt(T);
    
    // For the Pilot Simulations we need the : Sum of Y1_samples where Y1 := S1_T
    // Same for Y2_samples , where Y2 := S2_T as control variate
    double S1_T_sum = 0.0; 
    double S2_T_sum = 0.0;
    
    // Sum of X*Y1 from sampling in order to later estimate the Cov(X,Y1)
    double product1_sum = 0.0;
    // Same for the Cov(X,Y2)
    double product2_sum = 0.0;
    
    double disc_payoff_sum = 0.0;
    
    // Main loop for the path prices generation
    int i,j ;
    for( i=0; i<num_pilot; i++){
    	
    	// Generating each ε_i ~ N(0,1) for asset i 
    	// Here there are only 2 needed.
    	
    	double epsilon1 = Stand_Normal_Rand();
    	double epsilon2 = Stand_Normal_Rand();
    	
    	double S1_T = S1 * exp(nu_1 + v_1 * epsilon1);
    	double S2_T = S2 * exp(nu_2 + v_2 * epsilon2);
    	
    	//Setting up the payoff of the 2 assets Rainbow Price
    	double minimum1 = min(S1_T,S2_T);
    	
    	double disc_payoff = max(minimum1 - K,0.0) * exp( - r * T) ;
    	
    	// In general the formula of : Cov(X,Yi) = E[XY_k,i] - E[X]E[Y_k.i] with k=1,2 and i : 0 up to num_pilot
    	// So the 1st term is correspondingly per asset the :
    	double product1 = S1_T * disc_payoff ; 
    	double product2 = S2_T * disc_payoff;
    	// Accumulate ΣΥ1_i and ΣY2_i
    	S1_T_sum += S1_T;
    	S2_T_sum += S2_T;
    	// Accumulate Σ(X_i)
    	disc_payoff_sum += disc_payoff ;
    	// Accumulate Σ(X_i,Y1_i)
    	product1_sum += product1;
    	// Accumulate Σ(X_i,Y1_i)
    	product2_sum += product2;
    	
	}
	   // sample covariance estimator for Cov(X,Y) from pilot:
    // Cov(X,Y) ≈ [ Σ(XY) - (ΣX)(ΣY)/N ] / (N-1)  (PDF gives sample covariance formula for pilot) 
    double sample1_cov = (product1_sum - S1_T_sum * disc_payoff_sum / num_pilot) / (num_pilot - 1);
    double sample2_cov = (product2_sum - S2_T_sum * disc_payoff_sum / num_pilot) / (num_pilot - 1);

    // For this application the control variate is Yi = Si_T.
    // The PDF states E[S_T] = S0 e^{rT} and Var(S_T) = S0^2 e^{2rT}(e^{σ^2 T} - 1). 
    double VarY1 = S1 * S1 * exp(2.0 * r * T) * (exp(v1 * v1 * T) - 1.0); // Var(Yi) known in closed form 
    double VarY2 = S2 * S2 * exp(2.0 * r * T) * (exp(v2 * v2 * T) - 1.0); // Var(Y) known in closed form 
   
	
	double ExpY_1 = S1 * exp( r * T); // E[Y1] is known 
	double ExpY_2 = S2 * exp( r * T); // E[Y2] is known
	
	// Now we need to update the c* per asset where we obtained from minimizing Var(X_c)
	// X_c = X + c1 * ( Y1 - E[Y1]) + c2 * ( Y2 - E[Y2]) 
	// So : Var(X_c) = Var(X) + 2 * Σ ( c_i * Cov(X,Y_i) + ΣΣ c_i * c_j * Cov(Yi,Yj)
	
	// Setting up the c1* and c2*
	double c1 = -sample1_cov/VarY1;
	double c2 = -sample2_cov/VarY2;

	double control_var_sum = 0.0;
	double control_var_squared_sum = 0.0;
	
	// Now we proceed to main looping interations:
	for(int i=0; i<num_sims; i++){
		
		double epsilon1 = Stand_Normal_Rand();
		double epsilon2 = Stand_Normal_Rand();
		
		//Create the new S1_T,S2_T:
		double new_S1_T = S1 * exp( nu_1 + v_1*epsilon1 ); 
		double new_S2_T = S2 * exp( nu_2 + v_2*epsilon2 );
		
		double new_minimum = min( new_S1_T, new_S2_T);
		double new_disc_payoff = exp( -r * T) * max( new_minimum - K, 0.0);
		
		double control_var = new_disc_payoff + c1 * (new_S1_T - ExpY_1) + c2 * (new_S2_T - ExpY_2);
		
		control_var_sum += control_var;
		control_var_squared_sum += control_var * control_var;
	}
	
    double control_var_average = control_var_sum / num_sims;
    sample_var = (control_var_squared_sum - num_sims*control_var_average*control_var_average) / (num_sims - 1);
    
    return control_var_average;	
}

int main() {
	//unsigned int seed = 123456u; // pick any fixed integer seed
	//srand(seed);
	// Getting the actual parameters :
	
    int num_sims = 4000000;   // Number of simulated asset paths
    int num_pilot = 10000;
	double S1 = 100.0;      // Stock1 price
    double S2 = 100.0;      // Stock2 price
    double K = 100.0;      // Strike price
    double r = 0.05;       // Risk-free rate (5%)
    double v1 = 0.2;        // Volatility of the underlying (20%)
    double v2 = 0.2;        // Volatility of the underlying (20%)
    double T = 1.0;        // One year until expiration

    // Via the function we compute the Call price of Rainbow Option
    
    double call_price = MC_Call_Price_Rainbow_CV(S1,S2,K,r,v1,v2,T,num_sims,num_pilot);
    
    // Relative error for alpha = 5% (95% confidence)
    double z = inverse_of_normal_cdf(0.975, 0.0, 1.0);
    double rel_error = z * std::sqrt(sample_var / num_sims) / call_price;
    
    // Print out the results 
    cout << "Number of Paths for Control Variates: " << num_sims << endl;
    cout << "Number of Pilot sims:" << num_pilot << endl;
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