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
        double A1, A2, phi1, phi2, L1, L2, dispersion1, dispersion2,
                s1, s2, z_crit;
        std::vector<double> true_zs;

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
    = std::vector<std::string>{"A1", "A2", "phi1", "phi2", "L1", "L2",
                                "dispersion1", "dispersion2",
                                "s1", "s2", "z_crit"};

std::vector<double> MyModel::xs;
std::vector<double> MyModel::ys;
std::vector<double> MyModel::zs;
std::vector<double> MyModel::sig_zs;
std::vector<double> MyModel::vs;
std::vector<double> MyModel::sigmas;

inline MyModel::MyModel(RNG& rng)
:A1(800.0*rng.rand())
,A2(800.0*rng.rand())
,phi1(-M_PI + 2*M_PI*rng.rand())
,phi2(-M_PI + 2*M_PI*rng.rand())
,L1(2.0*rng.rand())
,L2(2.0*rng.rand())
,dispersion1(400.0*rng.rand())
,dispersion2(400.0*rng.rand())
,s1(-4.0 + 8.0*rng.rand())
,s2(-4.0 + 8.0*rng.rand())
,z_crit(-3.0 + 2.0*rng.rand())
,true_zs(sig_zs.size())
{
    for(size_t i=0; i<true_zs.size(); ++i)
    {
        true_zs[i] = zs[i] + sig_zs[i]*rng.randn();
    }
}

inline double MyModel::perturb(RNG& rng)
{
    int which = rng.rand_int(11);
    double logh = 0.0;

    if(which == 0)
    {
        A1 += 800.0*rng.randh();
        wrap(A1, 0.0, 800.0);
    }
    else if(which == 1)
    {
        A2 += 800.0*rng.randh();
        wrap(A2, 0.0, 800.0);
    }
    else if(which == 2)
    {
        phi1 += 2*M_PI*rng.randh();
        wrap(phi1, -M_PI, M_PI);
    }
    else if(which == 3)
    {
        phi2 += 2*M_PI*rng.randh();
        wrap(phi2, -M_PI, M_PI);
    }
    else if(which == 4)
    {
        L1 += 2.0*rng.randh();
        wrap(L1, 0.0, 2.0);
    }
    else if(which == 5)
    {
        L2 += 2.0*rng.randh();
        wrap(L2, 0.0, 2.0);
    }
    else if(which == 6)
    {
        dispersion1 += 400.0*rng.randh();
        wrap(dispersion1, 0.0, 400.0);
    }
    else if(which == 7)
    {
        dispersion2 += 400.0*rng.randh();
        wrap(dispersion2, 0.0, 400.0);
    }
    else if(which == 8)
    {
        s1 += 8.0*rng.randh();
        wrap(s1, -4.0, 4.0);
    }
    else if(which == 9)
    {
        s2 += 8.0*rng.randh();
        wrap(s2, -4.0, 4.0);
    }
    else if(which == 10)
    {
        z_crit += 2.0*rng.randh();
        wrap(z_crit, -3.0, -1.0);
    }
    else
    {
        int reps = 1 + rng.rand_int(10);
        for(int i=0; i<reps; ++i)
        {
            int k = rng.rand_int(true_zs.size());
            logh -= -0.5*pow((true_zs[k] - zs[k])/sig_zs[k], 2);
            true_zs[k] += sig_zs[k]*rng.randh();
            logh += -0.5*pow((true_zs[k] - zs[k])/sig_zs[k], 2);
        }
    }

    return logh;
}

inline double MyModel::log_likelihood() const
{
    double logl = 0.0;

    double A, phi, L, dispersion, s;
    double r, theta, mu, var;
    for(size_t i=0; i<xs.size(); ++i)
    {
        if(true_zs[i] < z_crit)
        {
            A = A1;
            phi = phi1;
            L = L1;
            dispersion = dispersion1;
            s = s1;
        }
        else
        {
            A = A2;
            phi = phi2;
            L = L2;
            dispersion = dispersion2;
            s = s2;
        }
        r = sqrt(xs[i]*xs[i] + ys[i]*ys[i]);

        // KV
        theta = atan2(ys[i], xs[i]);
        mu = A*sin(theta - phi);

        // KS
//        mu = A*(xs[i]*sin(phi) - ys[i]*cos(phi));

        // KF
//        mu = A*tanh((xs[i]*sin(phi) - ys[i]*cos(phi))/L);

        var = pow(dispersion*exp(s*r), 2) + pow(sigmas[i], 2);
        logl += -0.5*log(2*M_PI*var) - 0.5*pow(vs[i] - mu, 2)/var;
    }

    return logl;
}

inline std::vector<char> MyModel::to_blob() const
{
    std::vector<char> result(11*sizeof(double));
    auto pos = &result[0];
    for(double value: {A1, A2, phi1, phi2, L1, L2, dispersion1, dispersion2,
                        s1, s2, z_crit})
    {
        std::memcpy(pos, &value, sizeof(value));
        pos += sizeof(value);
    }
    return result;
}

inline void MyModel::from_blob(const std::vector<char>& blob)
{
    auto pos = &blob[0];
    for(double* value: {&A1, &A2, &phi1, &phi2, &L1, &L2, &dispersion1, &dispersion2, &s1, &s2, &z_crit})
    {
        std::memcpy(value, pos, sizeof(double));
        pos += sizeof(double);
    }
}


inline std::string MyModel::to_string() const
{
    return Tools::render(std::vector<double>{A1, A2, phi1, phi2, L1, L2, dispersion1, dispersion2, s1, s2, z_crit}, ",");
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

