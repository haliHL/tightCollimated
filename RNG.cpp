#include "RNG.hpp"

RNG::RNG ()
{        
        r_ = gsl_rng_alloc (gsl_rng_default);
        gsl_rng_set (r_, gsl_rng_default_seed);
}

RNG::RNG (unsigned long int seed)
{        
        r_ = gsl_rng_alloc (gsl_rng_default);
        gsl_rng_set (r_, seed);
}

RNG::~RNG ()
{
        gsl_rng_free (r_);
}

double RNG::get_uniform_rn (double xmin, double xmax) const
{
        return xmin + (xmax - xmin) * gsl_rng_uniform (r_); 
}

double RNG::get_gaussian_rn (double sigma) const
{
        return  gsl_ran_gaussian (r_, sigma); 
}

unsigned long int RNG::get_uniform_int () const
{
        return gsl_rng_get (r_);
}

unsigned long int RNG::get_poissonian_int ( double mean ) const
{
        return gsl_ran_poisson ( r_, mean );
}

unsigned long int RNG::get_binomial_int (double p, unsigned int n) const
{
	    return gsl_ran_binomial (r_, p, n);
}
