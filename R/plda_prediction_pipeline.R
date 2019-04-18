# assume rows represent samples
# n.samples <- nrow(dat)
# pct.training <- 0.7

# split data into training and testing
# idx.train <- sample(1:n.samples, round(n.samples*pct.training))
# data.train <- dat[idx.train, ]
# data.test <- dat[-idx.train, ]

# train pLDA on training
# fit.plda <- plda(dat, k=n.topics, lam)
# subset beta by cutoff
# retrain LDA on seleted genes, get beta* and gamma*
# use beta* to get gamma_test, train model on gamma*
# use trained model to predict on gamma_test