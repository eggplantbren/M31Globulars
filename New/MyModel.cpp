#include "MyModel.h"
#include "DNest4/code/DNest4.h"

Data MyModel::data("cher_data.csv");

MyModel::MyModel()
:A(num_components)
,phi(num_components)
,sigma(num_components)
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
    for(int i=0; i<num_components; ++i)
    {
        A[i] = 1000.0*rng.rand();
        phi[i] = 2.0*M_PI*rng.rand();
        sigma[i] = 1000.0*rng.rand();
    }
}

double MyModel::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    int which = rng.rand_int(3);
    int k = rng.rand_int(num_components);
    if(which == 0)
    {
        A[k] += 1000.0*rng.randh();
        DNest4::wrap(A[k], 0.0, 1000.0);
    }
    else if(which == 1)
    {
        phi[k] += 2.0*M_PI*rng.randh();
        DNest4::wrap(phi[k], 0.0, 2.0*M_PI);
    }
    else
    {
        sigma[k] += 1000.0*rng.randh();
        DNest4::wrap(sigma[k], 0.0, 1000.0);
    }

    return logH;
}

int MyModel::choose_component(double metallicity) const
{
    if(metallicity > -1.0)
        return 1;
    else
        return 0;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;
    for(size_t i=0; i<data.x.size(); ++i)
    {
        int component = choose_component(data.metallicity[i]);
        double mu = A[component]*sin(data.theta[i] - phi[component]);
        double var = pow(sigma[component], 2) + pow(data.verr[i], 2);

        logL += -0.5*log(2.0*M_PI*var) - 0.5*pow(data.v[i] - mu, 2)/var;
    }

    return logL;
}

void MyModel::print(std::ostream& out) const
{
    for(int i=0; i<num_components; ++i)
        out << A[i] << ' ' << phi[i] << ' ' << sigma[i] << ' ';
}

std::string MyModel::description() const
{
    return std::string("");
}

