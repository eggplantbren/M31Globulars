# Define the number of parameters
num_params = 3

# Define the names of the parameters.
# If you can't be bothered, use parameter_names=rep(NA, num_params).
parameter_names = c("A", "phi", "dispersion")

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

    params["A"]   = qunif(us[1], 0,  200)
    params["phi"] = qunif(us[2], -pi, pi)
    params["dispersion"] = qunif(us[3], 0, 400)

    return(params)
}

# Function that takes a vector of parameters and returns the
# log likelihood.
log_likelihood = function(params)
{
    mu = params["A"]*sin(data$theta - params["phi"])
    sig = sqrt(params["dispersion"]^2 + data$sig_v^2)
    sum(dnorm(data$v, mu, sig, log=TRUE))
}

