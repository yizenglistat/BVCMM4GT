preprocess_fake <- function(folder="./data", vars=NULL, linear=NULL, twoway=FALSE, N_test=600, reformat=FALSE)
{

	if(reformat){
		print("real raw data is protected due to privacy reasons! Hence formate always FALSE")
	}else{
		formated <- readRDS(glue("{folder}/formated_fake.RData"))
		S <- formated$S
		X <- formated$X
		Z <- formated$Z
		Y <- formated$Y

		u <- do_normalization(X[,1], log=TRUE)
		
		u_seq <- seq(-1.25, 1.25, length.out=N_test)
		X <- as.matrix(cbind(1, 2*X[,-1]-1))

		if(!is.null(vars)){
			X <- X[,vars]
		}

		if(twoway){
			X <- do_twoway(X)
		}

		if(!is.null(linear)){
			X_lin <- as.matrix(X[, linear])
			X <- X[, -linear]
		}else{
			X_lin <- NULL
		}
		
		N <- length(u)
		
		npools <- nrow(Z)
		pools <- Z[,2]

		# no varying intercept
		# X <- X[,-1]

		data <- list(testing="DT", X=as.matrix(X), u=u, X_lin=X_lin,
			N=N, npools=npools, pools=pools, N_test=N_test, 
			u_seq=u_seq, se=NULL, sp=NULL, S=S, nsites=length(unique(S)), Z=Z, Y=Y)

		return(data)
	}	
}






























