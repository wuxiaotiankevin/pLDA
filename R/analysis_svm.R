# input
# dat: raw data,
# y: label for raw data (cell types)
# props: proportion of training
# n_runs: number of svm repetition

# output
# svm predicted accuracy on current setting

# run once
run.SVM.once <- function(dat, y, train_prop, run_name, seed) {
  set.seed(seed)
  y <- factor(y)

  # split train and test
  tr_idx <- sort(sample(1:length(y), round(train_prop*length(y))))
  te_idx <- 1:length(y)
  te_idx <- te_idx[!(1:length(y) %in% tr_idx)]

  xtr <- dat[tr_idx, ]
  ytr <- y[tr_idx]

  xte <- dat[te_idx, ]
  yte <- y[te_idx]

  # run svm
  # print(ytr)
  # print(yte)
  model <- svm(x=xtr, y=ytr)
  pred <- predict(model, xte)

  # calculate accuracy
  accu <- sum(pred==yte)/length(yte)

  return(data.frame(run_name, train_prop, accu))
}

#' @export
runSVM <- function(dat, y, props, run_name, n_runs) {
  res <- data.frame()
  for (train_prop in props) {
    tmp <- mclapply(1:n_runs, function(x) run.SVM.once(dat, y, train_prop, run_name, x), mc.cores = n_runs)
    # res <- rbind(res, rbind.fill(tmp))
    res <- rbind(res, do.call("rbind", tmp))
  }
  res
}
