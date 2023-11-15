formate_fake <- function(folder='./data'){
	set.seed(1)
	data<-read.csv(glue("{folder}/fake.csv"),header=TRUE)

	##################################################################
	# Replace NA individual tests with 0 -- for data structuring only
	# NA individual test means pool tested negative, so contributing
	# members are diagnosed negative

	id = which(is.na(data$Ind.C.Result))
	data$Ind.C.Result[id] = 0



	##################################################################
	# Divide the data into Swab and Urine data sets
	# Note that Urine samples are tested individually (pool size = 1)

	pool.id<-data$Pool.ID
	data.gts<-data[!is.na(pool.id),]
	data.ind<-data[is.na(pool.id),]


	##################################################################
	# Creating the design matrices

	X.ind<-cbind(1,data.ind[,5:10])
	X.pool<-cbind(1,data.gts[,5:10])


	###################################################################
	# Building the Z and Y matrices

	Z.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=7)
	Y.ind.C<-matrix(-99,nrow=dim(data.ind)[1],ncol=4)

	Z.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=7)
	Y.pool.C<-matrix(-99,nrow=dim(data.gts)[1],ncol=4)


	##################################
	# For Swab testing

	id.ind<-1:(dim(data.gts)[1])

	pool.id<-unique(data.gts$Pool.ID)
	n<-length(pool.id)

	track.CT<-1

	for(i in 1:n){

	temp<-data.gts[data.gts$Pool.ID==pool.id[i],]
	temp.id<-id.ind[data.gts$Pool.ID==pool.id[i]]

	CT.res <- temp$Ind.C.Result

	cj<-length(CT.res)

	CT.retest<- sum(CT.res)>0

	Z.pool.C[track.CT,1]<-CT.retest
	Z.pool.C[track.CT,2]<-cj
	Z.pool.C[track.CT,3]<-1		## swab pool assay used
	Z.pool.C[track.CT,4:(cj+3)]<-temp.id

	Y.pool.C[temp.id,1]<-0
	Y.pool.C[temp.id,2]<-1
	Y.pool.C[temp.id,3]<-track.CT

	if(CT.retest==0){
	track.CT<-track.CT+1
	}

	if(CT.retest>0){
	tid<-(track.CT+1):(track.CT+cj)
	Z.pool.C[tid,1]<-CT.res
	Z.pool.C[tid,2]<-1
	Z.pool.C[tid,3]<-2		## swab individual assay used
	Z.pool.C[tid,4]<-temp.id

	Y.pool.C[temp.id,1]<-CT.res	## Diagnosed status
	Y.pool.C[temp.id,2]<-2
	Y.pool.C[temp.id,4]<-tid
	track.CT<-max(tid)+1
	}
	}


	Z.pool.C<-Z.pool.C[1:(track.CT-1),]

	ZnC<-dim(Z.pool.C)[1]
	YnC<-dim(Y.pool.C)[1]

	#################################
	# For Urine testing

	Z.ind.C[,1]<-data.ind$Ind.C.Result 	## Test outcome for Chlamydia
	Z.ind.C[,2]<-1                          ## Pool size    
	Z.ind.C[,3]<-3				## urine assay used
	Z.ind.C[,4]<-(1:(dim(data.ind)[1]))+YnC

	Y.ind.C[,1]<-data.ind$Ind.C.Result	## Diagnosed status
	Y.ind.C[,2]<-1
	Y.ind.C[,3]<-1:(dim(data.ind)[1])+ZnC


	###################################################
	# Putting everything together 

	X<-rbind(X.pool,X.ind)
	Z.C<-rbind(Z.pool.C,Z.ind.C)
	Y.C<-rbind(Y.pool.C,Y.ind.C)



	###################################################
	# Adding clinic ID to Y matrices
	Y.C = cbind(Y.C[,1:2], -99, Y.C[,3:4])

	Y.C[1:dim(data.gts)[1],3] = data.gts$Clinic.ID
	Y.C[(dim(data.gts)[1]+1):dim(Y.C)[1],3] = data.ind$Clinic.ID

	#####################################################################
	# Reenumerate Clinic ID -- needed for coding structure of routines

	clinic.id = unique(Y.C[,3])
	K = length(clinic.id)

	Y.C.new = Y.C

	for(i in 1:K){
		tmp.id = which(Y.C[,3] == clinic.id[i])
		Y.C.new[tmp.id,3] = i
	}

	Y.C = Y.C.new

	saveRDS(list(S=Y.C[,3], X=X[,-1], Y=Y.C[,-3], Z=Z.C), glue('{folder}/formated_fake.RData'))

}






