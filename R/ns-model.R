# Define the number of parameters
num_params = 7

# Define the names of the parameters.
# If you can't be bothered, use parameter_names=rep(NA, num_params).
parameter_names = c("A1", "phi1", "dispersion1",
                    "A2", "phi2", "dispersion2",
                    "z_crit")

# A dataset
# This could easily be loaded from an external file
data = read.table("../DNest5/data.txt")
data = data.frame(x=data[, 1], y=data[, 2],
                  z=data[, 3], sig_z = data[, 4],
                  v=data[, 5], sig_v = data[, 6])
data$theta = atan2(data$y, data$x)

# Function that takes a vector of Uniform(0, 1) variables
# and returns a vector of the actual parameters. This
# function implicitly defines the prior for the parameters.
us_to_params = function(us)
{
    # Vector to be returned as the result of the function
    params = rep(NA, num_params)

    # Apply the names
    names(params) = parameter_names

    #### You'll only need to edit this function below this line ####

    params["A1"]   = qunif(us[1], 0,  800)
    params["phi1"] = qunif(us[2], -pi, pi)
    params["dispersion1"] = qunif(us[3], 0, 400)
    params["A2"]   = qunif(us[4], 0,  800)
    params["phi2"] = qunif(us[5], -pi, pi)
    params["dispersion2"] = qunif(us[6], 0, 400)
    params["z_crit"] = qunif(us[7], -3, -1)

    return(params)
}

# Function that takes a vector of parameters and returns the
# log likelihood.
log_likelihood = function(params)
{
    mu = params["A1"]*sin(data$theta - params["phi1"])
    greater = data$z > params["z_crit"]
    mu[greater] = params["A2"]*sin(data$theta[greater] - params["phi2"])

    sig = sqrt(params["dispersion1"]^2 + data$sig_v^2)
    sig[greater] = sqrt(params["dispersion2"]^2 + data$sig_v[greater]^2)
    sum(dnorm(data$v, mu, sig, log=TRUE))
}

