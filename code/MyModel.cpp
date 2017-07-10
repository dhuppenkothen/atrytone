#include "MyModel.h"
#include "DNest4/code/Utils.h"
#include "Data.h"
#include "Lookup.h"
#include <cmath>
#include <stdexcept>
#include <gsl/gsl_sf_gamma.h>

using namespace std;
using namespace DNest4;

const Data& MyModel::data = Data::get_instance();
const int& nlines = Data::get_instance().get_nlines();
//const std::vector<double> line_pos = Data::get_instance().get_line_pos();
const PHAData& pha = Data::get_instance().get_pha();

// Initialise the static distribution
const DNest4::Cauchy MyModel::cauchy(0.0, 1.0);

MyModel::MyModel()
:dopplershift(2*nlines+1, 1, false, MyConditionalPrior())
,noise_normals(pha.bin_lo.size())
,mu(pha.counts.size())
{
}

double MyModel::gaussian_cdf(double x, double x0, double gamma)
{
	return 0.5*(1. + Lookup::erf((x-x0)/(gamma*sqrt(2.))));
}

template <typename ConstIntType, typename ConstFloatType,
            typename FloatType, typename IndexType, typename UIndexType>

void MyModel::rmf_fold(IndexType len_source, const ConstFloatType *source,
                IndexType len_num_groups, const ConstIntType *num_groups,
                IndexType len_first_chan, const ConstIntType *first_chan,
                IndexType len_num_chans, const ConstIntType *num_chans,
                IndexType len_response, const ConstFloatType *resp,
                IndexType len_counts, FloatType *counts,
                UIndexType offset)
{
    //int flag = 0;
//    if ( ( len_num_groups != len_source ) ||
//         ( len_first_chan != len_num_chans ))
//		throw RMFConvolutionFailure();

    // How many groups are in the current energy bin?
    IndexType current_num_groups = 0;

    // How many channels are in the current group?
    IndexType current_num_chans = 0;

    // What is the current channel of the output (counts) array?
    //register IndexType current_chan = 0;

    FloatType source_bin_ii;
    const ConstFloatType *resp_tmp = resp;
    const ConstIntType *first_chan_tmp = first_chan;
    const ConstIntType *num_chans_tmp = num_chans;
    FloatType *counts_tmp = counts;


    for (size_t ii = 0; ii < len_source; ii++ ) {

      // ii is the current energy bin
      source_bin_ii = source[ ii ];

      current_num_groups = num_groups[ ii ];

      while( current_num_groups ) {

//        if ( ( IndexType(first_chan_tmp - first_chan) >= len_num_chans ) ||
//             ( UIndexType(*first_chan_tmp) < offset ) )
//                //throw RMFConvolutionFailure();


        counts_tmp = counts + *first_chan_tmp - offset;
        current_num_chans = *num_chans_tmp;
        first_chan_tmp++;
        num_chans_tmp++;

//        if ( ( (IndexType(counts_tmp-counts) + current_num_chans) > len_counts )
//             ||
//             ( (IndexType(resp_tmp-resp) + current_num_chans) > len_response ) )
//                throw RMFConvolutionFailure();

        while ( current_num_chans ) {

          *counts_tmp += *resp_tmp * source_bin_ii;
          counts_tmp++;
          resp_tmp++;
          current_num_chans--;

        }
        current_num_groups--;

      }

    } // end for ii

  }



void MyModel::calculate_mu()
{ 

	// define line positions
	const vector<double>& line_pos = data.get_line_pos();

	// declare shifted line positions
	vector<double> line_pos_shifted(line_pos.size());

	const vector<double>& f_left = pha.bin_lo;
	const vector<double>& f_right = pha.bin_hi;
	const vector<double>& f_mid = pha.bin_mid;

	//cout<<"f_mid: "<<f_mid<<endl;

	// assign constant background to model
	mu.assign(mu.size(), 0.0);
	mu_bkg.assign(mu.size(), 0.0);

	// get amplitudes and widths from the RJObject 
	const vector< vector<double> >& dopplershiftcomponents = dopplershift.get_components();
 
	vector<double> amplitude(nlines), sign(nlines), width(nlines);
	double dshift;

//        // I'm only interested in a specific region of the spectrum
//        // right now, so let's only look at that!
	const double& f_min = data.get_f_min();
	const double& f_max = data.get_f_max();

//	// extra lines!
//	vector<double> e_line_pos;
//	e_line_pos.assign(3, 0.0);
//	e_line_pos[0] = 0.28325;
//	e_line_pos[1] = 0.28625;
//	e_line_pos[2] = 0.5270;
//
//
//        vector<double> e_line_sign;
//        e_line_sign.assign(3, 0.0);
//        e_line_sign[0] = -1. ;
//        e_line_sign[1] = 1.;
//        e_line_sign[2] = -1.;
//
        for(size_t j=0; j<dopplershiftcomponents.size(); j++)
        {

			dshift = dopplershiftcomponents[j][0];
	
			for (int i=0; i<nlines; i++)
				{
					line_pos_shifted[i] = line_pos[i] * (1./(1. + dshift));
					amplitude[i] = exp(dopplershiftcomponents[j][i+1]);
					//logq[i] = dopplershiftcomponents[j][i+1+nlines];
//      					sign[i] = dopplershiftcomponents[j][i+1+2*nlines];		
					width[i] = dopplershiftcomponents[j][i+1+nlines];
				}
	
        	        for(size_t i=0; i<mu.size(); i++)
        	        {
		                if (f_left[i] < f_min)                        
		                       continue;
		                if (f_right[i] > f_max)
	        	               continue;
	       		        else
	                	{       

					for (int k=0; k<nlines; k++)
						{
						// Integral over the Lorentzian distribution
//						if (sign[k] < dopplershift.get_conditional_prior().get_pp()) 
//							s = -1;
//						else 
//							s = 1;
						int s = -1;
						if ((std::abs(f_right[i] - line_pos_shifted[k]) < 5.*width[k]) && 
						   (std::abs(f_left[i] - line_pos_shifted[k]) < 5.*width[k])) 
							mu[i] += s*amplitude[k]*(gaussian_cdf(f_right[i], line_pos_shifted[k], width[k])
										- gaussian_cdf(f_left[i], line_pos_shifted[k], width[k]));
						}
                		}	 
			}

	}

        // fold through the ARF
        // code taken from sherpa
        for (size_t ii = 0; ii < mu.size(); ii++ )
		{ 
			// power law background
                        mu_bkg[ ii ] = exp(log(background/(slope + 1.)) + log(f_right[ ii ]) * (slope + 1.))
                                                 - exp(log(background/(slope + 1.)) + log(f_left[ ii ]) * (slope + 1.));


//			for (size_t jj = 0; jj< e_line_pos.size(); jj++)
//				{
//			// extra interstellar line
//			// this is a fudge until I can actually include proper physics models
//                        	mu_bkg[ ii ] += e_line_sign[ jj ] * e_amp[ jj ]*(gaussian_cdf(f_right[ ii ], e_line_pos[ jj ], e_width[ jj ])
//                                	            - gaussian_cdf(f_left[ ii ], e_line_pos[ jj ], e_width[ jj ]));
//                      		}
			// add Doppler-shifted lines to background, exponentiate and multiply by 
			// spectral response
			mu[ ii ] = exp(mu[ ii ] + log(mu_bkg[ ii ])) * pha.arf.specresp[ ii ];

		}


	vector<double> y(mu.size());
        double alpha = exp(-1./noise_L);

	// noise process could come both from the source or the detector!
 	// which is why I put it in between the ARF and the RMF
	for(size_t i=0; i<mu.size(); i++)
	{
                if(i == 0)
                        y[i] = noise_sigma*noise_normals[i];
                else
                        y[i] = alpha*y[i-1] + noise_sigma*noise_normals[i];

              //if ((f_left[i] < f_min) & (i > 0))
              //        y[i-1]=0.0;
              //else if (f_right[i] > f_max)
              //        y[i]=0.0;

                mu[i] *= exp(y[i]);

	}

        counts.assign(mu.size(), 0.0);

        rmf_fold(mu.size(), &mu[0], 
	  pha.rmf.n_grp.size(), &pha.rmf.n_grp[0],
	  pha.rmf.f_chan.size(), &pha.rmf.f_chan[0],
	  pha.rmf.n_chan.size(), &pha.rmf.n_chan[0],
	  pha.rmf.matrix.size(), &pha.rmf.matrix[0],
	  counts.size(), &counts[0],
	  pha.rmf.offset);

	counts.resize(f_left.size());
}

void MyModel::from_prior(RNG& rng)
{
	do
	{
		background = cauchy.generate(rng);
	}while(std::abs(background) > 25.0);
	background = exp(background);
	
//        e_amp.assign(3, 0.0);

//	for (size_t j=0; j < e_amp.size(); j++)
//	{
//		do
//		{
//	    	e_amp[j] = cauchy.generate(rng);
//		}while(std::abs(e_amp[j]) > 25.0);
//		e_amp[j] = exp(e_amp[j]);
//	}
//
//        e_width.assign(3, 0.0);
//
//        for (size_t j=0; j < e_width.size(); j++)
//        {
//                do
//                {
//                e_width[j] = cauchy.generate(rng);
//                }while(std::abs(e_width[j]) > 25.0);
//                e_width[j] = exp(e_width[j]);
//        }
//
//
	slope = -2.0 + 4.0 * rng.rand();
	
	dopplershift.from_prior(rng);

	noise_sigma = exp(log(1E-3) + log(1E3)*rng.rand());
	noise_L = exp(log(0.01*Data::get_instance().get_f_range())
                        + log(1000)*rng.rand());

	calculate_mu();

}

double MyModel::perturb(RNG& rng)
{
	double logH = 0.;

	int which;

	// A 50-50 decision to perturb a 'multivariate' thing or a 'univariate' thing
	if(rng.rand() <= 0.5)
	{
		which = rng.rand_int(2);

		if(which == 0)	// Proposal for the RJObject describing the lines
			logH += dopplershift.perturb(rng);
		else			// Proposal for normals for the OU process
		{
			if(rng.rand() <= 0.5)	// Propose to move only one
			{
				int i = rng.rand_int(noise_normals.size());
				logH -= -0.5*pow(noise_normals[i], 2);
				noise_normals[i] += rng.randh();
				logH += -0.5*pow(noise_normals[i], 2);



			}
			else					// Regenerate many
			{
                int reps = (int)pow(noise_normals.size(), rng.rand());
                for(int i=0; i<reps; ++i)
                {
                    int k = rng.rand_int(noise_normals.size());
                    noise_normals[k] = rng.randn();
                }

			}
		}
	}
	else
	{
		which = rng.rand_int(4);
		if(which == 0)
		{
            		background = log(background);
            		logH += cauchy.perturb(background, rng);
            		if(std::abs(background) > 25.0)
                		logH = -1E300;
            		background = exp(background);
		}

//                if(which == 1)
//                {
//			for (size_t j=0; j < e_amp.size(); j++)
//			{
//                        	e_amp[j] = log(e_amp[j]);
//                        	logH += cauchy.perturb(e_amp[j], rng);
//                       	 	if(std::abs(e_amp[j]) > 25.0)
//                         	       logH = -1E300;
//                        	e_amp[j] = exp(e_amp[j]);
//				}	
//
//                }
//
//                if(which == 2)
//                {
//			for (size_t j=0; j<e_width.size(); j++)
//			{
//                        	e_width[j] = log(e_width[j]);
//                        	logH += cauchy.perturb(e_width[j], rng);
//                        	if(std::abs(e_width[j]) > 25.0)
//                                	logH = -1E300;
//                       	 	e_width[j] = exp(e_width[j]);
//                	}
//		}
//

                else if(which == 1)
                        {
                        slope += 6.*rng.randh();
                        wrap(slope, -3., 3.);
                        }

		else if(which == 2)
		{
                        noise_sigma = log(noise_sigma);
                        logH += cauchy.perturb(noise_sigma, rng);
                        if(noise_sigma < -20.0 || noise_sigma > 0.0)
                                return -1E300;
                        noise_sigma = exp(noise_sigma);

		}
		else
		{
			noise_L = log(noise_L);
			noise_L += log(1E5)*rng.randh();
			wrap(noise_L, log(0.01*Data::get_instance().get_f_range()), log(1000.*Data::get_instance().get_f_range()));
			noise_L = exp(noise_L);
		}
	}

	calculate_mu();
	return logH;
}

double MyModel::log_likelihood() const
{

	const PHAData& pha = data.get_pha();

	const vector<double>& y1 = pha.counts;
	const vector<double>& f_left = pha.bin_lo;
        const vector<double>& f_right = pha.bin_hi;

	// I'm only interested in a specific region of the spectrum
	// right now, so let's only look at that!

        const double& f_min = data.get_f_min();
        const double& f_max = data.get_f_max();

        double logl = 0.;
	    for(size_t i=0; i<y1.size(); i++)
		{
                        if (f_left[i] < f_min)
                                continue;
                        if (f_right[i] > f_max)
                                continue;
			else
				{
					logl += -counts[i] + y1[i]*log(counts[i]) - gsl_sf_lngamma(y1[i] + 1.);
				}
 		}
	return logl;

}

void MyModel::print(std::ostream& out) const
{
        out.precision(25);

        const vector<double>& f_left = pha.bin_lo;
        const vector<double>& f_right = pha.bin_hi;

        // I'm only interested in a specific region of the spectrum
        // right now, so let's only look at that!

	const double& f_min = data.get_f_min();
        const double& f_max = data.get_f_max();

        out<<background<<' '<<slope<<' ';
//	for (size_t j=0; j< e_amp.size(); j++)
//		out<<e_amp[j]<<' ';
//
//        for (size_t j=0; j< e_amp.size(); j++)
//		out<<e_width[j]<<' ';


	out<<noise_L<<' '<<noise_sigma<<' ';
        
	dopplershift.print(out);

	for(size_t i=0; i<mu.size(); i++)
		{
                if (f_left[i] < f_min)
                                continue;
                if (f_right[i] > f_max)
                                continue;
                else
			out<<mu[i]<<' ';
 		}

        for(size_t i=0; i<counts.size(); i++)
                {
                if (f_left[i] < f_min)
                                continue;
                if (f_right[i] > f_max)
                                continue;
                else
                        out<<counts[i]<<' ';
                }

}

string MyModel::description() const
{
	return string("");
}

