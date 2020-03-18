#####Co-metric analysis
setwd("~/Google Drive/Manuscritos/MS_review_Through Dimensions of diversity/data_analysis")

####importing data####
metric_data<-read.table("matrixI_25-03-2019.txt",header=TRUE) #metrics by articles matrix 
metric_data<- metric_data[,-which(colSums(metric_data)==0)]

####co-occurence metric table####
co.metric<-as.matrix(metric_data)%*%t(as.matrix(metric_data)) #co-occurence metric
diag(co.metric)<-0 #setting the occurence of the metric with itselfe as zero

####Equivalence metric - calculated according to Callon et al. apud He#######
Ci_j.matrix<-matrix(nrow=60,ncol=60,byrow=TRUE, dimnames= list(rownames(co.metric), colnames(co.metric))) #number of co-occurence of pairs of metrics Cij  element in equation 3 of He and KD
for(i in 1:nrow(co.metric)) {       #numerator component of equation 3 in He and KD 
  Ci_j.matrix[i,]<-(co.metric[i,]^2)
}

metric.tot<-colSums(t(metric_data)) #total de vezes que cada m?trica ocorreu no conjunto de artigos
CixCj.matrix<-matrix(nrow=60,ncol=60,byrow=TRUE, dimnames= list(rownames(co.metric), colnames(co.metric))) #matrix to receive denominator of Eq.3 in He and Kd. Frequencies of each metric multiplied
for(i in 1:length(metric.tot)){		
	CixCj.matrix[i,]<-metric.tot[i]*metric.tot #denominator of Eq 3 in He and KD
}							
diag(CixCj.matrix)<-0 #filling the diagonal with zeros
equivalence.matrix<-matrix(NA,nrow=60,ncol=60,byrow=TRUE, dimnames= list(rownames(co.metric), colnames(co.metric))) #matrix to receive the equivalence values
for(i in 1:nrow(equivalence.matrix)){    
	equivalence.matrix[i,]<-((as.vector(Ci_j.matrix[i,]))/(as.vector(CixCj.matrix[i,]))) #calculating equivalence matrix according to Eq.3 in He and KD
}
equivalence.matrix #equivalence matrix
equival.matrix<-as.dist(equivalence.matrix) #tranforming in dist object 
equival.matrix

#####calculation of metric dendrogram########
cluster_metrics<- hclust(equival.matrix,method="ward") #direct calculated with Equivalence matrix, not recommended
quartz()
plot(cluster_metrics,labels=rownames(co.metric))
matrix.euclid<-dist(equivalence.matrix,method="euclidean") #transforming equivalence to square euclidean matrix as suggested by neff and corley
matrix.sqrtEuclidean<-matrix.euclid^2
cluster_metrics.SqrtEuclid<-hclust(matrix.sqrtEuclidean,method="ward.D") #cluster with square euclidean equivalence matrix using Ward agglomerative method
quartz()
par(mar=c(3,1,1,5)) 
plot(as.dendrogram(cluster_metrics.SqrtEuclid),horiz= T) #dendrogram with square euclidean and Ward method 

