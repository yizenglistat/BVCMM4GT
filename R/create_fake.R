#############################################################################
#############################################################################
##### This script will create a fake chlamydia data set to mimic the Iowa
##### group testing data set.
#############################################################################
#############################################################################
if(FALSE){
## Load necessary package ##
	library(mvtnorm)

	J = 65		## number of sites
	c = 4		## pool size
	N = 5000	## number of individuals

	## To perform Dorfman group testing algorithm ##
	Dorfman.decode.diff.error<-function(Y.true,Se,Sp,cj){
		N<-length(Y.true)
		Jmax<-N+N/cj
		J<-1

		Y<-matrix(-99,nrow=N, ncol=4)
		Z<-matrix(-99,nrow=Jmax,ncol=cj+3)


		for(j in 1:(N/cj)){
			prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se[1],1-Sp[1])
			Z[J,1]<-rbinom(n=1,size=1,prob=prob)
			Z[J,2]<-cj
			Z[J,3]<-1	## swab pool used
			Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
			Y[((j-1)*cj+1):(j*cj),2]<-1
			Y[((j-1)*cj+1):(j*cj),3]<-J
			J<-J+1
			if(Z[J-1,1]==1){
				for(k in ((j-1)*cj+1):(j*cj)){
					prob<-ifelse(Y.true[k]>0,Se[2],1-Sp[2])
					Z[J,1]<- rbinom(n=1,size=1,prob=prob)
					Z[J,2]<-1
					Z[J,3]<-2	## swab ind used
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


	## True parameter values ##
	Se.true = c(0.91,0.99)
	Sp.true = c(0.99,0.98)
	beta.true = lambda.true = rep(-99,7)
	beta.true[1] = -2
	beta.true[2] = -0.037
	beta.true[3] = -0.164
	beta.true[4] = 0.145
	beta.true[5] = 0.137
	beta.true[6] = 0.732
	beta.true[7] = 0.029
	lambda.true[1] = -0.25
	lambda.true[2] = 0.000
	lambda.true[3] = 0.000
	lambda.true[4] = 0.003
	lambda.true[5:7] = 0.000
	V.true = diag(lambda.true)
	A.true = diag(7)
	A.true[lower.tri(A.true)] = rnorm(7*6/2,0,0.5)
	D.true = V.true %*% A.true %*% t(A.true) %*% V.true
	b.true = t(rmvnorm(J, rep(0, 7), diag(7)))

	## Generate covariates ##
	set.seed(1)
	age = round(rnorm(N, 35, 7), 2)
	set.seed(NULL)
	race = rbinom(N, 1, 0.5)
	new_partner = rbinom(N, 1, 0.5)
	multi_partners = rbinom(N, 1, 0.13)
	contact = rbinom(N, 1, 0.05)
	symptoms = rbinom(N, 1, 0.28)
	X = cbind(age,race,new_partner,multi_partners,contact,symptoms)
	X = t((t(X) - apply(X,2,mean)) / apply(X,2,sd))
	X = cbind(1,X)


	## Generate individuals' true statuses ##
	Y.true = prob.save = rep(-99, N)
	sid = sample(1:J, N, replace = TRUE)	## site id for each individual
	for(i in 1:N){
			exb = exp(X[i,] %*% beta.true + X[i,] %*% diag(lambda.true) %*% A.true %*% b.true[,sid[i]])
	                prob = exb / (1 + exb)
			prob.save[i] = prob
	                Y.true[i] = rbinom(1, 1, prob)
	}

	## Perform Dorfman group testing ##
	group_data =  Dorfman.decode.diff.error(Y.true,Se.true,Sp.true,c)
	G = group_data$Z
	Y = group_data$Y
	Y = cbind(Y[,1:2],sid,Y[,3:4])	## adds site id
	Y[,1] = Y.true		## initializes as diagnose status
	Y1 = Y.true

	## Structure data similar to data application ##
	data.gt = data.frame()
	for(i in 1:N){
		site = Y[i,3]
		if(-99 %in% Y[i,4:5]){
			pid = Y[i,4]
			P.C.result = G[pid,1]
			if(P.C.result == 1){print("Error");break}
			if(P.C.result == 0){ind.res = NA}	## Does not complete individual test
			tmp = data.frame(i,site,pid,"Swab",age[i],race[i],new_partner[i],
					multi_partners[i],contact[i],symptoms[i],ind.res,0)
			names(tmp) = c("Patient","Clinic.ID","Pool.ID","Specimen","Age",
					"Race","New.Partner","Multiple.Partners","Contact","Symptoms",
					"Ind.C.Result","Pool.C.Result")
		} else if(!(-99 %in% Y[i,4:5])){
			pid = Y[i,4]
			P.C.result = G[pid,1]
			if(P.C.result == 0){print("Error");break}
			if(P.C.result == 1){
				pid2 = Y[i,5]
				ind.res = G[pid2,1]
			}
			tmp = data.frame(i,site,pid,"Swab",age[i],race[i],new_partner[i],
					multi_partners[i],contact[i],symptoms[i],ind.res,1)
			names(tmp) = c("Patient","Clinic.ID","Pool.ID","Specimen","Age",
					"Race","New.Partner","Multiple.Partners","Contact","Symptoms",
					"Ind.C.Result","Pool.C.Result")
		} else{
			print("Error")
			break
		}
		data.gt = rbind(data.gt,tmp)
		print(i)
	}


	## To add in Urine testing; i.e., individual level testing with pool size = 1 ##
	J = 65		## number of sites
	c = 1		## pool size
	N = 1000	## number of individuals

	Pool.test<-function(Y.true,Se,Sp,cj){
	N<-length(Y.true)
	Jmax<-N/cj 
	J<-1

	Y<-matrix(-99,nrow=N, ncol=3) 
	Z<-matrix(-99,nrow=Jmax,ncol=cj+3) 


	for(j in 1:(N/cj)){
	prob<-ifelse(sum(Y.true[((j-1)*cj+1):(j*cj)])>0,Se,1-Sp)
	Z[J,1]<-rbinom(n=1,size=1,prob=prob)
	Z[J,2]<-cj
	Z[J,3]<-3	## urine used
	Z[J,4:(cj+3)]<-((j-1)*cj+1):(j*cj)
	Y[((j-1)*cj+1):(j*cj),2]<-1
	Y[((j-1)*cj+1):(j*cj),3]<-J
	J<-J+1
	}

	J<-J-1
	Z<-Z[1:J,]

	return(list("Z"=Z, "Y"=Y))
	}


	## True parameter values ##
	Se.true = c(0.98)
	Sp.true = c(0.99)

	## Generate covariates ##
	set.seed(1)
	age = round(rnorm(N, 35, 7), 2)
	set.seed(NULL)
	race = rbinom(N, 1, 0.5)
	new_partner = rbinom(N, 1, 0.5)
	multi_partners = rbinom(N, 1, 0.13)
	contact = rbinom(N, 1, 0.05)
	symptoms = rbinom(N, 1, 0.28)
	X = cbind(age,race,new_partner,multi_partners,contact,symptoms)
	X = t((t(X) - apply(X,2,mean)) / apply(X,2,sd))
	X = cbind(1,X)

	## Generate individuals' true statuses ##
	Y.true = rep(-99, N)
	sid = sample(1:J, N, replace = TRUE)	## site id for each individual
	for(i in 1:N){
			exb = exp(X[i,] %*% beta.true + X[i,] %*% diag(lambda.true) %*% A.true %*% b.true[,sid[i]])
	                prob = exb / (1 + exb)
	                Y.true[i] = rbinom(1, 1, prob)
	}

	## Perform testing ##
	group_data =  Pool.test(Y.true,Se.true,Sp.true,c)
	G = group_data$Z
	Y = group_data$Y
	Y = cbind(Y[,1:2],sid,Y[,3])	## adds site id
	Y[,1] = Y.true		## initializes as diagnose status

	## Structure data similar to data application ##
	data.ind = data.frame()
	for(i in 1:N){
		site = Y[i,3]
		pid = Y[i,4]
		ind.res = G[pid,1]
		tmp = data.frame(i,site,"NA","Urine",age[i],race[i],new_partner[i],
				multi_partners[i],contact[i],symptoms[i],ind.res,NA)	## NA pool, no pooled test
		names(tmp) = c("Patient","Clinic.ID","Pool.ID","Specimen","Age",
				"Race","New.Partner","Multiple.Partners","Contact","Symptoms",
				"Ind.C.Result","Pool.C.Result")
		data.ind = rbind(data.ind,tmp)
		print(i)
	}

	## Combine data and save as excel file ##
	data = rbind(data.gt,data.ind)
	write.csv(data, file = "data/fake.csv", row.names=FALSE)
}


