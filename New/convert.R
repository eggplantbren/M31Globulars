# Convert the output to Cher's preferred CSV format.

posterior_samples = read.table("posterior_sample.txt", header=FALSE, skip=1)

colnames(posterior_samples) = c("A2", "phi2", "sigma2",
                                "A1", "phi1", "sigma1",
                                "M_crit")

df = data.frame(sigma1=posterior_samples$sigma1,
                A1=posterior_samples$A1,
                phi1=posterior_samples$phi1,
                sigma2=posterior_samples$sigma2,
                A2=posterior_samples$A2,
                phi2=posterior_samples$phi2,
                M_crit=posterior_samples$M_crit)

rownames(df) = NULL
write.csv(df, "postdf.csv", row.names=FALSE)
