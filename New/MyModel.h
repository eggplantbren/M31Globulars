#ifndef DNest4_Template_MyModel
#define DNest4_Template_MyModel

#include "DNest4/code/DNest4.h"
#include <ostream>
#include <vector>
#include "Data.h"

class MyModel
{
    private:

        // Dataset
        static Data data;

        // Number of components
        static constexpr int num_components = 2;

        // Amplitude, orientation angle, and velocity dispersion per component
        std::vector<double> A, phi, sigma;

        // Metallicity cut
        double m_crit;

        // Standard normals for the true metallicities
        std::vector<double> ns;

        // Choose the appropriate component
        int choose_component(double true_metallicity,
                             double reported_metallicity) const;


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
};

#endif

