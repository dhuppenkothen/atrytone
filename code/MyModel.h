#ifndef _MyModel_
#define _MyModel_

#include "DNest4/code/RNG.h"
#include "MyConditionalPrior.h"
#include <ostream>
#include "Data.h"
#include "DNest4/code/Distributions/Cauchy.h"
#include "DNest4/code/RJObject/RJObject.h"
#include <vector>
#include <stdexcept>


 
class MyModel
{
	private:
		// Reference to the data
		static const Data& data;

        // A useful cauchy distribution
        static const DNest4::Cauchy cauchy;

		// A flat background level
		double background;

		// The Lorentzians
		DNest4::RJObject<MyConditionalPrior> dopplershift;
		// Extra white noise on the flux
		std::vector<double> noise_normals;

		double noise_sigma, noise_L;

		// Poisson mean
		// these should probably be long doubles
		// need to fix that some time!
		std::vector<double> mu;
                std::vector<double> counts;

//		std::vector<double> mu_small;

		// Calculate mu from scratch:
		void calculate_mu();


		template <typename ConstIntType, typename ConstFloatType,
            		typename FloatType, typename IndexType, typename UIndexType>

		// Fold through RMF:
		void rmf_fold(IndexType len_source, const ConstFloatType *source,
                IndexType len_num_groups, const ConstIntType *num_groups,
                IndexType len_first_chan, const ConstIntType *first_chan,
                IndexType len_num_chans, const ConstIntType *num_chans,
                IndexType len_response, const ConstFloatType *resp,
                IndexType len_counts, FloatType *counts,
                UIndexType offset);

 
	public: 
		// Constructor only gives size of params
		MyModel();

		// Generate the point from the prior
		void from_prior(DNest4::RNG& rng);

		// Metropolis-Hastings proposals
		double perturb(DNest4::RNG& rng);

		// Likelihood function
		double log_likelihood() const;

		// Print to stream
		void print(std::ostream& out) const;

		// Return string with column information
		std::string description() const;

		static double gaussian_cdf(double x, double x0, double gamma);
};

#endif

