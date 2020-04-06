#ifndef RNG_HPP
#define RNG_HPP

extern "C" {
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
}

class RNG {
        public:
                // The default constructor generates a random number
                // generator that is seeded by time.
                RNG ();
                RNG (unsigned long int seed);
                ~RNG ();

                // making these random number extraction functions const is not
                // strictly speaking correct because it changes the state of
                // the random number generator. It is statistically correct
                // though because the statistical properties of the next random
                // number are unchanged.
                double get_uniform_rn ( double xmin, double xmax ) const;
                double get_gaussian_rn ( double sigma ) const;
                unsigned long int get_uniform_int () const;
                unsigned long int get_poissonian_int ( double mean ) const;
                unsigned long int get_binomial_int (double p, unsigned int n) const;

        private :

                gsl_rng *r_;

};


#endif

