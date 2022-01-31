#ifndef DNest5_Template_MyModel_hpp
#define DNest5_Template_MyModel_hpp

#include <cmath>
#include <fstream>
#include <memory>
#include <UniformModel.hpp>
#include <Tools/Misc.hpp>

namespace DNest5_Template
{

using DNest5::ParameterNames, Tools::RNG, Tools::wrap;

class MyModel
{
    private:

        double sigma, A1, A2, phi1, phi2, L, s, M_crit;
        double errorbar_multiplier;
        std::vector<double> true_zs;

        // Data
        static std::vector<double> xs, ys, zs, sig_zs, vs, sigmas;

    public:
        MyModel(RNG& rng);
        void us_to_params();
        double log_likelihood() const;
        double perturb(RNG& rng);
        std::vector<char> to_blob() const;
        void from_blob(const std::vector<char>& blob);
        std::string to_string() const;
        static const ParameterNames parameter_names;
        static void load_data();

};

/* Implementations follow */

std::vector<double> MyModel::xs;
std::vector<double> MyModel::ys;
std::vector<double> MyModel::zs;
std::vector<double> MyModel::sig_zs;
std::vector<double> MyModel::vs;
std::vector<double> MyModel::sigmas;

const ParameterNames MyModel::parameter_names
    = std::vector<std::string>{"sigma", "A1", "A2", "phi1", "phi2", "L",
                               "s", "M_crit", "errorbar_multiplier"};

MyModel::MyModel(RNG& rng)
{
    sigma = pow(10.0, 2.0 + 2.0*rng.randn());
    A1 = sigma*pow(10.0, 2.0*rng.randn());
    A2 = sigma*pow(10.0, 2.0*rng.randn());
    phi1 = 2.0*M_PI*rng.rand();
    phi2 = 2.0*M_PI*rng.rand();
    L = pow(10.0, 1.0 + 2.0*rng.randn());
    s = pow(10.0, 2.0*rng.randn());
    M_crit = -3.0 + 3.5*rng.rand();
    errorbar_multiplier = rng.randn(); // If < 0.5 it's 1, otherwise pareto
    for(size_t i=0; i<zs.size(); ++i)
        true_zs.emplace_back(zs[i] + sig_zs[i]*rng.randn());
}

double MyModel::perturb(RNG& rng)
{
    double logh = 0.0;
    int which = rng.rand_int(10);
    if(which == 0)
    {
        sigma = log10(sigma);
        logh -= -0.5*pow((sigma - 2.0)/2.0, 2);
        sigma += 2.0*rng.randh();
        logh += -0.5*pow((sigma - 2.0)/2.0, 2);
        sigma = pow(10.0, sigma);
    }
    else if(which == 1)
    {
        A1 = log10(A1/sigma);
        logh -= -0.5*pow(A1/2.0, 2);
        A1 += 2.0*rng.randh();
        logh += -0.5*pow(A1/2.0, 2);
        A1 = sigma*pow(10.0, A1);
    }
    else if(which == 2)
    {
        A2 = log10(A2/sigma);
        logh -= -0.5*pow(A2/2.0, 2);
        A2 += 2.0*rng.randh();
        logh += -0.5*pow(A2/2.0, 2);
        A2 = sigma*pow(10.0, A2);
    }
    else if(which == 3)
    {
        phi1 += 2.0*M_PI*rng.randh();
        wrap(phi1, 0.0, 2.0*M_PI);
    }
    else if(which == 4)
    {
        phi2 += 2.0*M_PI*rng.randh();
        wrap(phi2, 0.0, 2.0*M_PI);
    }
    else if(which == 5)
    {
        L = log10(L);
        logh -= -0.5*pow((L - 1.0)/2.0, 2);
        L += 2.0*rng.randh();
        logh += -0.5*pow((L - 1.0)/2.0, 2);
        L = pow(10.0, L);
    }
    else if(which == 6)
    {
        s = log10(s);
        logh -= -0.5*pow(s/2.0, 2);
        s += 2.0*rng.randh();
        logh += -0.5*pow(s/2.0, 2);
        s = pow(10.0, s);
    }
    else if(which == 7)
    {
        M_crit += 3.5*rng.randh();
        wrap(M_crit, -3.0, 0.5);
    }
    else if(which == 8)
    {
        errorbar_multiplier += rng.randh();
        wrap(errorbar_multiplier);
    }
    else if(which == 9)
    {
        int k = rng.rand_int(true_zs.size());
        logh -= -0.5*pow((true_zs[k] - zs[k])/sig_zs[k], 2);
        true_zs[k] += sig_zs[k]*rng.randh();
        logh += -0.5*pow((true_zs[k] - zs[k])/sig_zs[k], 2);
    }

    return logh;
}

double MyModel::log_likelihood() const
{
    double logl = 0.0;

    double A, phi, var, r, mu;
//  variance=(params["sigma"]*exp(-(newdata$r/(params["s"]*params["L"]))))^2+(newdata$V_sig)^2
//  line = rep(0, n)
//  subset = m > params["M_crit"]
//  line[subset] = params["A1"]*tanh((newdata$x[subset]*sin(params["phi1"])-newdata$y[subset]*cos(params["phi1"]))/params["L"])
//  line[!subset] = params["A2"]*tanh((newdata$x[!subset]*sin(params["phi2"])-newdata$y[!subset]*cos(params["phi2"]))/params["L"])
//  logL = sum(dnorm(newdata$V_M31, line, sd=sqrt(variance), log=TRUE))
//  

    // Mixture of delta(1) and pareto(xmin=1, alpha=1)
    double multiplier = 1.0;
    if(errorbar_multiplier > 0.5)
    {
        double u = 2.0*(errorbar_multiplier - 0.5);
        double e = -log(1.0 - u);
        multiplier = exp(e);
    }

    for(size_t i=0; i<xs.size(); ++i)
    {
        r = sqrt(pow(xs[i], 2) + pow(ys[i], 2));
        var = pow(sigma*exp(-r/(s*L)), 2) + pow(sigmas[i]*multiplier, 2);
        if(true_zs[i] < M_crit)
        {
            A = A1;
            phi = phi1;
        }
        else
        {
            A = A2;
            phi = phi2;
        }

        mu = A*tanh((xs[i]*sin(phi) - ys[i]*cos(phi))/L);
        logl += -0.5*log(2*M_PI*var) - 0.5*pow(vs[i] - mu, 2)/var;
    }

    return logl;
}

std::vector<char> MyModel::to_blob() const
{
    std::stringstream ss;
    ss.write(reinterpret_cast<const char*>(&sigma), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&A1), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&A2), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&phi1), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&phi2), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&L), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&s), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&M_crit), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&errorbar_multiplier), sizeof(double));

    std::string s = ss.str();
    std::vector<char> result(s.size());
    std::memcpy(&result[0], &s[0], s.size());

    return result;
}

void MyModel::load_data()
{
    std::fstream fin("data.txt", std::ios::in);
    double x, y, m, sig_z, v, sig;
    xs.clear();
    ys.clear();
    zs.clear();
    sig_zs.clear();
    vs.clear();
    sigmas.clear();

    while(fin >> x && fin >> y && fin >> m && fin >> sig_z && fin >> v && fin >> sig)
    {
        xs.push_back(x);
        ys.push_back(y);
        zs.push_back(m);
        sig_zs.push_back(sig_z);
        vs.push_back(v);
        sigmas.push_back(sig);
    }
    std::cout << "Loaded " << xs.size() << " data points." << std::endl;


    fin.close();
}

void MyModel::from_blob(const std::vector<char>& blob)
{
    std::stringstream ss;
    for(char c: blob)
        ss << c;
    ss.read(reinterpret_cast<char*>(&sigma), sizeof(double));
    ss.read(reinterpret_cast<char*>(&A1), sizeof(double));
    ss.read(reinterpret_cast<char*>(&A2), sizeof(double));
    ss.read(reinterpret_cast<char*>(&phi1), sizeof(double));
    ss.read(reinterpret_cast<char*>(&phi2), sizeof(double));
    ss.read(reinterpret_cast<char*>(&L), sizeof(double));
    ss.read(reinterpret_cast<char*>(&s), sizeof(double));
    ss.read(reinterpret_cast<char*>(&M_crit), sizeof(double));
    ss.read(reinterpret_cast<char*>(&errorbar_multiplier), sizeof(double));
}

std::string MyModel::to_string() const
{
    std::stringstream ss;
    ss << sigma << ',';
    ss << A1 << ',';
    ss << A2 << ',';
    ss << phi1 << ',';
    ss << phi2 << ',';
    ss << L << ',';
    ss << s << ',';
    ss << M_crit << ',';
    ss << errorbar_multiplier;
    return ss.str();
}

} // namespace

#endif

