#ifndef DNest5_Template_MyModel_hpp
#define DNest5_Template_MyModel_hpp

#include <cmath>
#include <fstream>
#include <UniformModel.hpp>
#include <Tools/Misc.hpp>

namespace DNest5_Template
{

using DNest5::ParameterNames, Tools::RNG, Tools::wrap;

class MyModel
{
    private:
        double A, phi, dispersion;

        // Data
        static std::vector<double> xs, ys, zs, sig_zs, vs, sigmas;

    public:
        inline MyModel(RNG& rng);
        inline double perturb(RNG& rng);
        inline double log_likelihood() const;
        inline std::vector<char> to_blob() const;
        inline void from_blob(const std::vector<char>& blob);
        inline std::string to_string() const;

        static ParameterNames parameter_names;
        static void load_data();
};

/* Implementations follow */

ParameterNames MyModel::parameter_names
    = std::vector<std::string>{"A", "phi", "dispersion"};

std::vector<double> MyModel::xs;
std::vector<double> MyModel::ys;
std::vector<double> MyModel::zs;
std::vector<double> MyModel::sig_zs;
std::vector<double> MyModel::vs;
std::vector<double> MyModel::sigmas;

inline MyModel::MyModel(RNG& rng)
:A(1000.0*rng.rand())
,phi(2*M_PI*rng.rand())
,dispersion(1000.0*rng.rand())
{

}

inline double MyModel::perturb(RNG& rng)
{
    int which = rng.rand_int(3);

    if(which == 0)
    {
        A += 1000.0*rng.randh();
        wrap(A, 0.0, 1000.0);
    }
    else if(which == 1)
    {
        phi += 2*M_PI*rng.randh();
        wrap(phi, 0.0, 2*M_PI);
    }
    else
    {
        dispersion += 1000.0*rng.randh();
        wrap(dispersion, 0.0, 1000.0);
    }

    return 0.0;
}

inline double MyModel::log_likelihood() const
{
    double logl = 0.0;

    double theta, mu, var;
    for(size_t i=0; i<xs.size(); ++i)
    {
        theta = atan2(ys[i], xs[i]);
        mu = A*sin(theta - phi);
        var = pow(dispersion, 2) + pow(sigmas[i], 2);
        logl += -0.5*log(2*M_PI*var) - 0.5*pow(vs[i] - mu, 2)/var;
    }

    return logl;
}

inline std::vector<char> MyModel::to_blob() const
{
    std::vector<char> result(3*sizeof(double));
    auto pos = &result[0];
    for(double value: {A, phi, dispersion})
    {
        std::memcpy(pos, &value, sizeof(value));
        pos += sizeof(value);
    }
    return result;
}

inline void MyModel::from_blob(const std::vector<char>& blob)
{
    auto pos = &blob[0];
    for(double* value: {&A, &phi, &dispersion})
    {
        std::memcpy(value, pos, sizeof(double));
        pos += sizeof(double);
    }
}


inline std::string MyModel::to_string() const
{
    return Tools::render(std::vector<double>{A, phi, dispersion}, ",");
}


void MyModel::load_data()
{
    std::fstream fin("../DNest5/data.txt", std::ios::in);
    double x, y, z, sig_z, v, sig;
    xs.clear();
    ys.clear();
    zs.clear();
    sig_zs.clear();
    vs.clear();
    sigmas.clear();

    while(fin >> x && fin >> y && fin >> z && fin >> sig_z && fin >> v && fin >> sig)
    {
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(z);
        sig_zs.push_back(sig_z);
        vs.push_back(v);
        sigmas.push_back(sig);
    }
    std::cout << "Loaded " << xs.size() << " data points." << std::endl;


    fin.close();
}

} // namespace

#endif

