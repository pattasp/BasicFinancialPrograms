#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstdlib>

using namespace std;

//Problem Set 2 : Question 4 | i) StopLossPut 

// Generate a random number between 0 and 1
// Return a uniform number in [0,1].
double UnifRand()
{
    return rand() / double(RAND_MAX);
}

// The inverse of the cumulative distribution function of a Normal
// distribution with mean mu and standard deviation sigma, used to generate
// gaussian random numbers - necessary for the Monte Carlo method below
double inverse_of_normal_cdf(double p, double mu, double sigma)
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

    if (std::abs(q) <= .425)
    {
        r = .180625 - q * q;
        val =
            q * (((((((r * 2509.0809287301226727 +
                      33430.575583588128105) * r +
                     67265.770927008700853) * r +
                    45921.953931549871457) * r +
                   13731.693765509461125) * r +
                  1971.5909503065514427) * r +
                 133.14166789178437745) * r +
                3.387132872796366608) /
            (((((((r * 5226.495278852854561 +
                  28729.085735721942674) * r +
                 39307.89580009271061) * r +
                21213.794301586595867) * r +
               5394.1960214247511077) * r +
              687.1870074920579083) * r +
             42.313330701600911252) * r + 1);
    }
    else
    {
        if (q > 0) r = 1 - p;
        else       r = p;

        r = std::sqrt(-std::log(r));

        if (r <= 5)
        {
            r += -1.6;
            val =
                (((((((r * 7.7454501427834140764e-4 +
                      .0227238449892691845833) * r +
                     .24178072517745061177) * r +
                    1.27045825245236838258) * r +
                   3.64784832476320460504) * r +
                  5.7694972214606914055) * r +
                 4.6303378461565452959) * r +
                1.42343711074968357734) /
                (((((((r * 1.05075007164441684324e-9 +
                      5.475938084995344946e-4) * r +
                     .0151986665636164571966) * r +
                    .14810397642748007459) * r +
                   .68976733498510000455) * r +
                  1.6763848301838038494) * r +
                 2.05319162663775882187) * r + 1);
        }
        else
        {
            r += -5;
            val =
                (((((((r * 2.01033439929228813265e-7 +
                      2.71155556874348757815e-5) * r +
                     .0012426609473880784386) * r +
                    .026532189526576123093) * r +
                   .29656057182850489123) * r +
                  1.7848265399172913358) * r +
                 5.4637849111641143699) * r +
                6.6579046435011037772) /
                (((((((r * 2.04426310338993978564e-15 +
                      1.4215117583164458887e-7) * r +
                     1.8463183175100546818e-5) * r +
                    7.868691311456132591e-4) * r +
                   .0148753612908506148525) * r +
                  .13692988092273580531) * r +
                 .59983220655588793769) * r + 1);
        }

        if (q < 0.0) val = -val;
    }

    return mu + sigma * val;
}

// Generate a gaussian random number from N(0,1)
double Standard_Normal_Rand()
{
    double u = UnifRand();
    while (u <= 0 || u >= 1) u = UnifRand();
    return inverse_of_normal_cdf(u, 0, 1);
}

// Global variable for the Sample Variance 	
double sample_var = 0.0 ;

// Stop-Loss naive method function of Hedging a Put Option via Monte Carlo Simulation

// Kindly reminder : Put = - max(St-K,0.0).

double StopLossPut(double S,double K,double mu,double r,double T,double v,int N,int num_sims){
	
	//Again we redifine the needed parameters
	
	double dt = T/N ;
	double nudt = (mu- 0.5* v* v)*dt;
	double vdt = v * sqrt(dt);
	
	double disc_payoff_sum = 0.0 ;
	double disc_payoff_squared_sum = 0.0;
	
	// Setting up the lattice==array to keep track of each path for the hedging strategy :
	
	double SPath[N+1]; // * we could also use a double vector == 2D dynamic array to prevent 
	// memory overfilling for N-> large 
	int Covered; 
	// Here we start again the loops to set up each possible payoffs : Pathing Generation 
	
	int i;
	for(i=0; i<num_sims; i++){
		double disc_payoff = 0;
		SPath[0] = S;
		//Note : The Stop-Loss needs to keep track when the option is ITM and when OTM:
		// So we need another variable to keep track of that situtation. For shortance we
		// call it "Covered" 
		//int Covered; // note to me : Should I declare it outside of for loop ?
		
		// 2 possible states : ITM -> Buy 1 unit of stock "Covered Position"
		// OTM -> Sell that previously bought stock to hold a "Naked Position" opposite of Covered
		// So we can use the Covered variable all the time : when Covered == 1 then ITM, Covered==0 OTM
		
		if(SPath[0] < K){  // Here SPath[0] < K indicates we have a Put Option 
			Covered = 1;
			disc_payoff += SPath[0]; // note to mslf ; why not negative ? 
		}else {
			Covered = 0;
		}
		
		// Now we make another for loop to track all the possible time steps :
		int j;
		for(j=0; j<N ; j++){
			double epsilon = Standard_Normal_Rand();
			SPath[j+1] = SPath[j] * exp(nudt + vdt* epsilon);
			//Again we check the ITM,OTM for each possible node:
			if(Covered == 1 && SPath[j+1] >= K){
				Covered =0;
				disc_payoff -= exp(-r*(j+1)*dt) * SPath[j+1]; //note to msfl : again here why not possitive ? 
			}else if(Covered == 0 && SPath[j+1] < K){
				Covered = 1;
				disc_payoff += exp(-r*(j+1)*dt) * SPath[j+1];
				
			}
		}
		//Now we need to check if the last Node at N-th index the Option is Excersiced !
		if(SPath[N] < K){
			disc_payoff -= exp(-r*T) * K ;	// again here negative instead of positive sign 
		}
		double disc_payoff_squared = disc_payoff*disc_payoff;
		disc_payoff_sum += disc_payoff;
		disc_payoff_squared_sum += disc_payoff_squared;
	}
	double disc_payoff_average = - disc_payoff_sum / num_sims ; //We use the negative sign "-" because we need to pay
	sample_var = (disc_payoff_squared_sum - num_sims*disc_payoff_average*disc_payoff_average)/(num_sims -1);
	return disc_payoff_average; 
}

int main(){
	
	//Setting up the actual MC simulation variables 
	double S = 27;
	double K = 30;
	double mu = 0.2 ;
	double r = 0.05 ;
	double T = 0.5;
	double v = 0.3;
	int N = 180;
	int num_sims = 13000;

	//The cost is always the alternate sign of the MC so the cost is positive :
    double cost =  StopLossPut(S,K,mu,r,T,v,N,num_sims);
	
	cout << "The cost for the StopLoss Hedging is: " << cost << endl;
	return 0;
}