#include <iostream>
#include <cmath>
#include <string>
#include <sstream>
#include <cstdlib>
#include <algorithm>

using namespace std;

// Problem Set 2 | Question 4 ii)

// Our task is to construct a Delta Hedging Put function
// that computes the cost of delta 'Δ' hedging for a EU Put Option
// again the input variables will be the same as the StopLoss Put

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

// Because of the nature of Delta give by the Black Scholes, we know that 
// Δ = Φ(d1) for call option || Δ = Φ(d1) - 1 for put option

double Phi(double x){
	return std::erfc(-x/std::sqrt(2.0)) / 2.0;
}

// Global variable for the Sample Variance 	
double sample_var = 0.0 ;

double MC_DeltaHedgingPut(double S,double K, double mu,double r,double T,double v,int N,int num_sims){
	
	double dt = T/N ;
	double nudt = (mu- 0.5* v* v)*dt;
	double vdt = v * sqrt(dt);
	
	double disc_payoff_sum = 0.0 ;
	double disc_payoff_squared_sum = 0.0;
	
	double SPath[N+1];
	double Position ;
	double delta;
	double d1; //this is from the Black Scholes PDE where d1 = log(S/K) + (r-0.5*v^2)*T / v*sqrt(T)
	
	double disc_payoff;
	
	//Initializing the Path possible payoffs via the for loop
	for(int i=0; i<num_sims; i++){
		SPath[0] = S;
		
		//Here d1 s= log(S0/K) + (r-0.5*v^2)*T / v*sqrt(T) but S_0 is the initial price so 
		// it must be equal to SPath[0]
		
		d1 = (log(SPath[0]/K) + (r+0.5*v*v)*T ) / (v*sqrt(T));
		
		// delta below is : Δ^i_0 -1 because we are calculating a Put Position
		delta = Phi(d1) - 1; //How many shares I have to hold in my position in order to Hedge
		
		// Here we are at the 1st path 
		disc_payoff =  -delta * SPath[0];
		
		Position = delta; // Keeping track of how many delta's I need between time steps 
		
		
		for(int j=0; j<N-1; j++){
			double epsilon = Standard_Normal_Rand();
			SPath[j+1] = SPath[j] * exp(nudt + vdt* epsilon);
			
			// We need to always keep track the time progression up until T so this must 
			// be implemented in d1 and thus inside the Delta.
			d1 = (log(SPath[j+1]/K) + (r+0.5*v*v) * (T - (j+1)*dt) ) / ( v*sqrt(T-(j+1)*dt));
			
			// At different nodes we need to see if the number of delta stocks needs to be revaluated
			delta = Phi(d1) - 1.0;
			
			// Now we need to calculate the currents delta's present value cashflow and 
			// add them to the previous value.
			
			disc_payoff += exp(-r*(j+1)*dt)*(Position -delta)*SPath[j+1];
			
			// Update based on our new delta the also new position we will need to hold
			Position = delta;
		}
		
		// Last step to maturity
		double epsilon = Standard_Normal_Rand();
		SPath[N] = SPath[N-1] * exp(nudt + vdt*epsilon);
		
		if(SPath[N] < K){
			//The Deltas for any time are in the range of [0,1]. So they can be e.g 0.45 
			// But at the maturity T the option is either Excercised or Not 
			// So at the end when the counterparty is going to excercise it will ask that 
			// 1 share he previously agreed on the strike price K. So we need to match it in 
			// our Hedging Strategy 
			
			//Also at t=T the difference between the time_steps is 0. So for d1 to not be infinity 
			// We calculate that time's payoff by hand.
			
			// Last we set it as Position + 1 because for Put -> Counterparty will want to "sell" 
			// So we need to be able to "buy" 1 stock in order to Hedge
			disc_payoff += exp(-r*T) * ((Position + 1.0) * SPath[N] - K);
		}else{ 
		    // We are in the OTM occassion so we wont need to hold any number of shares
			disc_payoff += exp(-r*T)* Position* SPath[N];
		}
		
		double disc_payoff_squared = disc_payoff * disc_payoff ;
		disc_payoff_sum += disc_payoff;
		disc_payoff_squared_sum += disc_payoff_squared;
		
	}
	
	double disc_payoff_average = -disc_payoff_sum / num_sims ;
	sample_var  = (disc_payoff_squared_sum - num_sims*disc_payoff_average*disc_payoff_average)/(num_sims-1);
	return disc_payoff_average; 
}


int main(){
	//Setting up our variables :
	
	double S = 27.0;
	double K = 30.0;
	double mu = 0.2;
	double r  = 0.05;
	double T = 0.5;
	double v = 0.3;
	int N = 180;
	int num_sims = 130000;
	
	double cost = MC_DeltaHedgingPut(S,K,mu,r,T,v,N,num_sims);
	cout << "The cost for Delta Hedging a Put Option is: " << cost << endl;
	
	return 0;
}