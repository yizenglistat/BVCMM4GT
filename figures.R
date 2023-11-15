# ------------------------------------- figures
pool_sizes <- c(5, 10)
folder <- "output"
model_names <- c('m1', 'm2')
sigma <- 0.5
nreps <- 5
testings <- c("IT", "DT", "AT")
knowns <- c(TRUE, FALSE)
Ns <- c(3000, 5000)

for(known in knowns){
	for(N in Ns){
		for(model_name in model_names){
			df <- c()
			for(pool_size in pool_sizes){
				for(testing in testings){
					print(glue("processing posterior {testing} data ..."))
					res <- post_summary(folder=folder,
						known=ifelse(known, "known", "unknown"),
						sample_size=N, 
						model_name=model_name, 
						pool_size=pool_size, testing=testing, 
						sigma=sigma, nreps=nreps, 
						cache=TRUE, replace=FALSE)
					is_pts <- is.na(res$levs)
					knots <- res$levs[!is_pts]
					med <- res$med[!is_pts]
					true <- res$true[!is_pts]
					lower <- res$lower[!is_pts]
					upper <- res$upper[!is_pts]
					df_testing <- data.frame(
						knots=knots,
						med=med,
						true=true,
						lower=lower,
						upper=upper,
						label=rep(c("f0","f1","f2"),each=length(knots)/3),
						testing=testing,
						pool=pool_size
					)
					df <- rbind(df, df_testing)
					print(glue("{testing} data done!"))
				}
			}
			rdata_file <- glue("{folder}/{ifelse(known, 'known', 'unknown')}/{N}/{model_name}/df_pool.RData")
			saveRDS(df, rdata_file)
		}
	}
}
