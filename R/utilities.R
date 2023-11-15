misc_str_pad <- function(string) return(format(round(as.numeric(string),3), nsmall=3))

misc_str <- function(string){
      if(is.na(string)) return(NA)
      if(grepl("\\(", string)){
            slice_start <- gregexpr("\\(", string)[[1]][1] + 1
            slice_end <- nchar(string) - 1

            return(glue(misc_str_pad(substr(string, 1, slice_start-2)), 
                  "(",
                  misc_str_pad(substr(string, slice_start, slice_end)),
                  ")"))
      }else{
            return(misc_str_pad(string))
      }
}

glue_tab_item <- function(x, y, digits=3) return(unlist(lapply(glue("{round(x,digits)}({round(y,digits)})"), misc_str)))


do_logit <- function(x) log(x/(1-x))
do_logit_inv <- function(x) return(1/(1+exp(-x)))

do_normalization <- function(x, log=FALSE)
{
	if(log){
		out <- (log(x)-3)
	}else{
		out <- (x-mean(x))/sd(x)
	}

	return(out)
}

do_unnormalization <- function(x, x_mean=24.22094, x_sd=6.712391, log=FALSE)
{
	if(log){
		out <- exp(x+3)
	}else{
		out <- x_sd*x+x_mean
	}

	return(out)
}

do_slice_add <- function(vec, p)
{
	nvec <- length(vec)
	nknots <- nvec/p
	slices <- seq(1, (p+1)*nknots, by=nknots)
	slice_base <- slices[1]:(slices[2]-1)
	base <- vec[slice_base]
	out <- c(base)
	for(idx in 2:p){
		slice <- slices[idx]:(slices[idx+1]-1)
		out <- c(out, base + vec[slice])
	}
	return(out)
}

do_twoway <- function(X){
	p <- ncol(X)
	names <- colnames(X)
	for(i in 2:(p-1)){
		for(j in (i+1):p){
			names <- c(names, glue::glue("{names[i]}_{names[j]}"))
			X <- cbind(X, as.vector(X[,i]*X[,j]))
		}
	}
	colnames(X) <- names
	return(X)
}

replace_id <- function(old, new, vec){
	vec[vec %in% old] <- new[match(vec, old, nomatch = 0)]
	return(vec)
}

remove_na <- function(df, col)
{
	ids <- df[is.na(df[,col]),c("seq_id","pool_id")]
	remove_pools <- ids$pool_id[!is.na(ids$pool_id)]
	df <- df[!(df$pool_id %in% remove_pools | df$seq_id %in% ids$seq_id),]

	return(df)
}


model_set <- function(model_name="m1", lambda=1.0){

	model_name <- substr(model_name,1,2)

	if(model_name=="m1"){
		
		f0 <- function(x, offset=-3.3) sin(pi*x/3) + offset
		f1 <- function(x, offset=-0.5) x^3/8 + offset
		f2 <- function(x, offset=0) -x^2/4 + offset
		

		f_list <- list(f0=f0, f1=f1, f2=f2)
	}

	if(model_name=="m2"){
		f0 <- function(x, offset=-3.3) cos(pi*x/2) + offset
		f1 <- function(x, offset=0.5) -x^3/8 + offset
		f2 <- function(x, offset=-0.7) x^2/4 + offset

		f_list <- list(f0=f0, f1=f1,f2=f2)
	
	}

	# x=seq(-3,3,0.001)
	# par(mfrow=c(1,3), pty='s')
	# plot(x,f0(x),'l')
	# plot(x,f1(x),'l')
	# plot(x,f2(x),'l')


	return(f_list)

}

do_IT<-function(Y_true,Se,Sp,cj=1){
   N<-length(Y_true)
   Jmax<-N/cj 
   J<-1

   Y<-matrix(-99,nrow=N, ncol=3) 
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


   for(j in 1:(N/cj)){
	   prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
	   Z[J,1]<-rbinom(n=1,size=1,prob=prob)
	   Z[J,2]<-cj
	   Z[J,3]<-2
	   Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
	   Y[((j-1)*cj+1):(j*cj),2]<-1
	   Y[((j-1)*cj+1):(j*cj),3]<-J
	   J<-J+1
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

do_MPT<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N/cj 
   J<-1

   Y<-matrix(-99,nrow=N, ncol=3) 
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


   for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

do_DT<-function(Y_true,Se,Sp,cj){
   N<-length(Y_true)
   Jmax<-N+N/cj
   J<-1

   Y<-matrix(-99,nrow=N, ncol=4)
   Z<-matrix(-99,nrow=Jmax,ncol=cj+3)


   for(j in 1:(N/cj)){
      prob<-ifelse(sum(Y_true[((j-1)*cj+1):(j*cj)])>0,Se[1],1-Sp[1])
      Z[J,1]<-rbinom(n=1,size=1,prob=prob)
      Z[J,2]<-cj
      Z[J,3]<-1   ## swab pool used
      Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
      Y[((j-1)*cj+1):(j*cj),2]<-1
      Y[((j-1)*cj+1):(j*cj),3]<-J
      J<-J+1
      if(Z[J-1,1]==1){
         for(k in ((j-1)*cj+1):(j*cj)){
            prob<-ifelse(Y_true[k]>0,Se[2],1-Sp[2])
            Z[J,1]<- rbinom(n=1,size=1,prob=prob)
            Z[J,2]<-1
            Z[J,3]<-2   ## swab ind used
            Z[J,4]<-k
            Y[k,2]<-2
            Y[k,4]<-J
            J<-J+1
         }
      }
   }

   J<-J-1
   Z<-Z[1:J,]

   return(list("Z"=Z, "Y"=Y))
}

do_AT<-function(Y_true, Se, Sp, cj){

	N<-length(Y_true)
	Jmax<-2*N/cj +N
	J<-1
	AT<-N/(cj^2)

	Y<-matrix(-99,nrow=N, ncol=5) 
	Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 

	Y.A<-array(-99,c(cj,cj,AT))
	ID.A<-array(-99,c(cj,cj,AT))
	ind<-1
	for(i in 1:AT){
	for(j in 1:cj){
	for(m in 1:cj){
	Y.A[m,j,i]<-Y_true[ind]
	ID.A[m,j,i]<-ind
	ind<-ind+1
	}}}



	for(s in 1:AT){

	array.yk<-Y.A[,,s]
	array.id<-ID.A[,,s]

	a<-rep(0,nrow(array.yk))
	b<-rep(0,ncol(array.yk))

	for(i in 1:cj){
	   prob<- ifelse(sum(array.yk[i,])>0, Se[1], 1-Sp[1])
	   g<- rbinom(n=1,size=1,prob=prob)
	   a[i]<-g
	   Z[J,1]<-g 
	   Z[J,2]<-cj 
	   Z[J,3]<-1
	   Z[J,4:(cj+3)]<-array.id[i,]
	   Y[array.id[i,],2]<-2
	   Y[array.id[i,],3]<-J
	   J<-J+1
	}
	for(j in 1:cj){
	   prob<- ifelse(sum(array.yk[,j])>0, Se[1], 1-Sp[1])
	   g<- rbinom(n=1,size=1,prob=prob)
	   b[j]<-g
	   Z[J,1]<-g 
	   Z[J,2]<-cj 
	   Z[J,3]<-1
	   Z[J,4:(cj+3)]<-array.id[,j]
	   Y[array.id[,j],4]<-J
	   J<-J+1
	}


	if(sum(a)>0 & sum(b)>0){
	array.yk1<-as.matrix(array.yk[(a==1),(b==1)])
	array.id1<-as.matrix(array.id[(a==1),(b==1)])
	for(i in 1:nrow(array.yk1)){
	for(j in 1:ncol(array.yk1)){
	   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
	   g<- rbinom(n=1,size=1,prob=prob)
	   Z[J,1]<-g 
	   Z[J,2]<-1 
	   Z[J,3]<-2
	   Z[J,4]<-array.id1[i,j]
	   Y[array.id1[i,j],2]<-3
	   Y[array.id1[i,j],5]<-J
	   J<-J+1
	}}}



	if(sum(a)>0 & sum(b)==0){
	array.yk1<-as.matrix(array.yk[(a==1),])
	array.id1<-as.matrix(array.id[(a==1),])
	for(i in 1:nrow(array.yk1)){
	for(j in 1:ncol(array.yk1)){
	   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
	   g<- rbinom(n=1,size=1,prob=prob)
	   Z[J,1]<-g 
	   Z[J,2]<-1 
	   Z[J,3]<-2
	   Z[J,4]<-array.id1[i,j]
	   Y[array.id1[i,j],2]<-3
	   Y[array.id1[i,j],5]<-J
	   J<-J+1
	}}}

	if(sum(a)==0 & sum(b)>0){
	array.yk1<-as.matrix(array.yk[,(b==1)])
	array.id1<-as.matrix(array.id[,(b==1)])
	for(i in 1:nrow(array.yk1)){
	for(j in 1:ncol(array.yk1)){
	   prob<- ifelse(array.yk1[i,j]>0, Se[2], 1-Sp[2])
	   g<- rbinom(n=1,size=1,prob=prob)
	   Z[J,1]<-g 
	   Z[J,2]<-1 
	   Z[J,3]<-2
	   Z[J,4]<-array.id1[i,j]
	   Y[array.id1[i,j],2]<-3
	   Y[array.id1[i,j],5]<-J
	   J<-J+1
	}}}

	} 

	J<-J-1
	Z<-Z[1:J,]

	return(list("Z"=Z, "Y"=Y))
}

do_stats <- function(est, std, true){

	est[is.na(est)] <- 0
	std[is.na(std)] <- 0


	bias <- apply(est-true, 1, mean)
	rmse <- apply(est-true, 1, function(err) rmse(err,0))
	ssd <- apply(est, 1, sd)
	ese <- apply(std, 1, mean)

	med <- apply(est, 1, quantile, probs=c(0.5))
	lower <- apply(est, 1, quantile, probs=c(0.025))
	upper <- apply(est, 1, quantile, probs=c(1-0.025))

	return(list(bias=bias, rmse=rmse, ese=ese, ssd=ssd, med=med, lower=lower, upper=upper))
}

post_summary <- function(folder='./output', sample_size=5000, model_name='m1', pool_size=5, testing='DT', sigma=0.5, nreps=500, cache=TRUE, replace=FALSE)
{
	cache_name <- glue('{folder}/N{sample_size}/{model_name}/cj{pool_size}/{testing}_summary.RData')
	if(cache&(!replace)&file.exists(cache_name)){
		return(readRDS(cache_name))
	}else{

		file_name <- glue('{folder}/N{sample_size}/{model_name}/cj{pool_size}/{testing}_summary.h5')
		
		if(file.exists(file_name)) file.remove(file_name)
		
		param <- H5File$new(file_name, mode = "a")		
		
		# create store matrix
		param_med <- param$create_dataset(name="param_med", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))
		param_std <- param$create_dataset(name="param_std", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))
		param_ext <- param$create_dataset(name="param_ext", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))
		param_eci <- param$create_dataset(name="param_eci", dtype=h5types$H5T_NATIVE_FLOAT, 
			space=H5S$new(dims=c(1, nreps), maxdims=c(Inf, nreps)),chunk_dims=c(1, nreps))

		SE <- c(0.95, 0.98)
		SP <- c(0.98, 0.99)

		pb <- progress::progress_bar$new(format="[:bar] :current/:total[:percent] elapsed :elapsed eta :eta",
	                       total = nreps)

		for(rep in 1:nreps){
			
			rdata_file <- glue("{folder}/N{sample_size}/{model_name}/cj{pool_size}/{testing}/fitted_task{rep}.RData")
			fitted <- readRDS(rdata_file)

			# sigma
			param_med[1,rep] <- fitted$sigma_hat_med
			param_std[1,rep] <- fitted$sigma_hat_std
			param_eci[1,rep] <- 1.0*I(fitted$sigma_hat_lower<sigma && sigma<fitted$sigma_hat_upper)

			# SE1
			param_med[2,rep] <- fitted$se_hat_med[1]
			param_std[2,rep] <- fitted$se_hat_std[1]	
			param_eci[2,rep] <- 1.0*I(fitted$se_hat_lower<SE[1] && SE[1]<fitted$se_hat_upper[1])

			# SE2
			param_med[3,rep] <- fitted$se_hat_med[2]
			param_std[3,rep] <- fitted$se_hat_std[2]	
			param_eci[3,rep] <- 1.0*I(fitted$se_hat_lower<SE[2] && SE[2]<fitted$se_hat_upper[2])

			# SP1
			param_med[4,rep] <- fitted$sp_hat_med[1]
			param_std[4,rep] <- fitted$sp_hat_std[1]	
			param_eci[4,rep] <- 1.0*I(fitted$sp_hat_lower<SP[1] && SP[1]<fitted$sp_hat_upper[1])

			# SP2
			param_med[5,rep] <- fitted$sp_hat_med[2]
			param_std[5,rep] <- fitted$sp_hat_std[2]	
			param_eci[5,rep] <- 1.0*I(fitted$sp_hat_lower<SP[2] && SP[2]<fitted$sp_hat_upper[2])


			if(model_name %in% c("m1", "m2")){
				f_list <- model_set(model_name=model_name)
				true <- c(sigma, SE, SP)
			}

			# beta(s) & f(s)
			binary <- 0
			N_test <- fitted$N_test
			nvar <- 3
			fitted$levs_test_hat_med <- c(seq(-3,3,length.out=N_test), seq(-3,3,length.out=N_test), seq(-3,3,length.out=N_test))
			binary_id <- 0
			
			for(var_id in 1:nvar){
				
				slice_start <- binary_id + 1 + N_test*(var_id-binary_id-1)
				slice_end <- binary_id + N_test*(var_id-binary_id)
				slice <- slice_start:slice_end
				
				true <- c(true, f_list[[glue("f{var_id-binary-1}")]](fitted$levs_test_hat_med[slice]))

				param_med[slice+5,rep] <- fitted$beta_test_hat_med[slice]
				param_std[slice+5,rep] <- fitted$beta_test_hat_std[slice]
				param_ext[var_id,rep] <- rmse(fitted$beta_test_hat_med[slice], true[slice])

			}

			param_ext[nvar+1,rep] <- fitted$time
			param_ext[nvar+2,rep] <- fitted$ntests_hat_med
			
			pb$tick()

		}

		results <- do_stats(param_med[,], param_std[,], true)

		results$levs <- c(NA, rep(NA, 4), fitted$levs_test_hat_med)

		results$true <- true
		results$param_med <- param_med[,]
		results$param_std <- param_std[,]
		results$param_eci <- param_eci[,]
		results$param_ext <- param_ext[,]

		param$close_all()
		
		saveRDS(results, cache_name)

		return(results)
	}
}

rmse <- function(y_hat, y_true) sqrt(mean((y_hat-y_true)^2))
















