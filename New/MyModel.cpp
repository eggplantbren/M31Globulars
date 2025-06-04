#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <iomanip>

Data MyModel::data("cher_data_fe_h_err.csv");

MyModel::MyModel()
:A(num_components)
,phi(num_components)
,sigma(num_components)
,L(num_components)
,ns(data.x.size())
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
    for(int i=0; i<num_components; ++i)
    {
        A[i] = 1000.0*rng.rand();
        phi[i] = -M_PI + 2.0*M_PI*rng.rand();
        sigma[i] = 1000.0*rng.rand();
        L[i] = exp(log(1E-1) + log(1E3)*rng.rand());
    }

    m_crit = -2.8 + 2.3*rng.rand();

    for(double& n: ns)
        n = rng.randn();
}

double MyModel::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(6);
    int k = rng.rand_int(num_components);
    if(which == 0)
    {
        A[k] += 1000.0*rng.randh();
        DNest4::wrap(A[k], 0.0, 1000.0);
    }
    else if(which == 1)
    {
        phi[k] += 2.0*M_PI*rng.randh();
        DNest4::wrap(phi[k], -M_PI, M_PI);
    }
    else if(which == 2)
    {
        sigma[k] += 1000.0*rng.randh();
        DNest4::wrap(sigma[k], 0.0, 1000.0);
    }
    else if(which == 3)
    {
        L[k] = log(L[k]);
        L[k] += log(1E3)*rng.randh();
        DNest4::wrap(L[k], log(1E-1), log(100.0));
        L[k] = exp(L[k]);
    }
    else if(which == 4)
    {
        m_crit += 2.3*rng.randh();
        DNest4::wrap(m_crit, -2.8, -0.5);
    }
    else
    {
        logH += DNest4::perturb_ns(ns, rng);
    }

    return logH;
}

int MyModel::choose_component(double true_metallicity,
                              double reported_metallicity) const
{
    // Model 1
    //    return 0;

    // Model 2.1
    if(reported_metallicity < 1000)
    {
        if(true_metallicity < m_crit)
            return 0;
        else
            return 1;
    }
    else if(std::abs(reported_metallicity - 1001) <= 1E-6)
        return 0;
    else
        return 1;


//    // Model 2.2
//    if(reported_metallicity < 1000)
//    {
//        if(true_metallicity < m_crit)
//            return 0;
//        else
//            return 1;
//    }
//    else if(std::abs(reported_metallicity - 1001) <= 1E-6)
//        return 1;
//    else
//        return 0;

    // Model 3
//    if(reported_metallicity < 1000)
//    {
//        if(true_metallicity < m_crit)
//            return 0;
//        else
//            return 1;
//    }
//    else if(std::abs(reported_metallicity - 1001) <= 1E-6)
//        return 0;
//    else
//        return 2;


}

double MyModel::log_likelihood() const
{
    double logL = 0.0;
    for(size_t i=0; i<data.x.size(); ++i)
    {
        double true_metallicity = data.metallicity[i]
                                    + data.metallicity_err[i]*ns[i];
        int component = choose_component(true_metallicity,
                                         data.metallicity[i]);

        // V model
        double mu = A[component]*sin(data.theta[i] - phi[component]);

        // S model
//        double mu = A[component]*(data.x[i]*sin(phi[component]) -
//                                  data.y[i]*cos(phi[component]));

        // F model
//        double mu = A[component]*tanh((data.x[i]*sin(phi[component]) -
//                                       data.y[i]*cos(phi[component]))/L[component]);

        double var = pow(sigma[component], 2) + pow(data.verr[i], 2);

        logL += -0.5*log(2.0*M_PI*var) - 0.5*pow(data.v[i] - mu, 2)/var;
    }

    return logL;
}

void MyModel::print(std::ostream& out) const
{
    out << std::setprecision(10);
    for(int i=0; i<num_components; ++i)
        out << A[i] << ' ' << phi[i] << ' ' << sigma[i] << ' ';
    out << m_crit << ' ';
}

std::string MyModel::description() const
{
    return std::string("");
}

