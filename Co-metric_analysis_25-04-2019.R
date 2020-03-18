#####Co-metric analysis
setwd("~/Google Drive/Manuscritos/MS_review_Through Dimensions of diversity/data_analysis")

####importing data####
metric_data<-read.table("matrixI_25-04-2019.txt", header=TRUE)
####co-occurence metric table####
co.metric<-as.matrix(metric_data)%*%t(as.matrix(metric_data)) #co-occurence metric
diag(co.metric)<-0 #setting the occurence of the metric with itselfe as zero

####Equivalence metric - calculated according to Callon et al. apud He#######
Ci_j.matrix<-matrix(nrow= length(rownames(co.metric)),ncol= length(rownames(co.metric)),byrow=TRUE, dimnames= list(rownames(co.metric), colnames(co.metric))) #number of co-occurence of pairs of metrics Cij  element in equation 3 of He and KD
for(i in 1:nrow(co.metric)) {       #numerator component of equation 3 in He and KD 
  Ci_j.matrix[i,]<-(co.metric[i,]^2)
}

metric.tot<-colSums(t(metric_data)) #total de vezes que cada m?trica ocorreu no conjunto de artigos
CixCj.matrix<-matrix(nrow=length(rownames(co.metric)),ncol=length(rownames(co.metric)),byrow=TRUE, dimnames= list(rownames(co.metric), colnames(co.metric))) #matrix to receive denominator of Eq.3 in He and Kd. Frequencies of each metric multiplied
for(i in 1:length(metric.tot)){		
	CixCj.matrix[i,]<-metric.tot[i]*metric.tot #denominator of Eq 3 in He and KD
}							
diag(CixCj.matrix)<-0 #filling the diagonal with zeros
equivalence.matrix<-matrix(NA,nrow=length(rownames(co.metric)),ncol=length(rownames(co.metric)),byrow=TRUE, dimnames= list(rownames(co.metric), colnames(co.metric))) #matrix to receive the equivalence values
for(i in 1:nrow(equivalence.matrix)){    
	equivalence.matrix[i,]<-((as.vector(Ci_j.matrix[i,]))/(as.vector(CixCj.matrix[i,]))) #calculating equivalence matrix according to Eq.3 in He and KD
}
equivalence.matrix #equivalence matrix
equival.matrix<-as.dist(equivalence.matrix) #tranforming in dist object 

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
plot(as.dendrogram(cluster_metrics.SqrtEuclid))
plot(ape::as.phylo(cluster_metrics.SqrtEuclid))

#######group density########

equival.matrix #matrix de equival??ncia
match(c("FEve", "Fdis","Fdiv", "FRic"),rownames(as.matrix(equival.matrix)))
#FEve, Fric, FDis, FEsp, FDiv

dens.struct<-(sum(equival.matrix[5,8],equival.matrix[5,9],equival.matrix[5,3],equival.matrix[5,19],
                  matrix.resu[8,9],matrix.resu[8,3],matrix.resu[8,19],matrix.resu[9,3],
                  matrix.resu[9,19],matrix.resu[3,19]))/(factorial(5)/(factorial(2)*factorial(3))) 
#calculo de densidade para o grupo filogen??tico
dens.phylo<-sum(matrix.resu[43,44],matrix.resu[43,41],matrix.resu[43,42],matrix.resu[44,41],
                matrix.resu[44,42],matrix.resu[41,42])/(factorial(4)/(factorial(2)*factorial(2)))
#calculo densidade para grupo funcional
dens.func<-(sum(matrix.resu[14,15],matrix.resu[14,12],
                matrix.resu[14,11],matrix.resu[14,13],matrix.resu[15,12],matrix.resu[15,11],
                matrix.resu[15,13],matrix.resu[12,11],matrix.resu[12,13],matrix.resu[11,13]))/
  (factorial(5)/(factorial(2)*factorial(3)))
#calculo de densidade para grupo total div, taxonomic entropy e total complex.
match(c("total.diversity", "taxonomic.entropy.phylogenetic"),rownames(as.matrix(equival.matrix)))
dens.total<-(sum(as.matrix(equival.matrix)[26,25]))/
  (factorial(2)/(factorial(2)*factorial(1)))
dens.total
#dominance e Hilsenhoff
dens.dominance<-(sum(matrix.resu[26,27]))/(factorial(2)/(factorial(2)*factorial(0)))
dens.dominance
#Rao e CWM
dens.rao<-(sum(matrix.resu[28,60]))/(factorial(2)/(factorial(2)*factorial(0)))
#tax. diversity e indice funcional n???o especificado
dens.taxdiversity<-
  dens.taxdiversity
#variance in tax. diversity e average tax distinctness
dens.variancetax<-(sum(matrix.resu[23,24]))/(factorial(2)/(factorial(2)*factorial(0)))
dens.variancetax
#berger park, rarity e menhinick
dens.berger<-(sum(matrix.resu[38,7],matrix.resu[38,39],matrix.resu[7,39]))/
  (factorial(3)/(factorial(2)*factorial(1)))
dens.berger
#margalef, phylogenetic distinctness, rarefaction e N1 diversity
dens.margalef<-(sum(matrix.resu[34,35],matrix.resu[34,33],matrix.resu[34,49],
                    matrix.resu[35,33],matrix.resu[35,49],matrix.resu[33,49]))/
  (factorial(4)/(factorial(2)*factorial(2)))
dens.margalef
#convex volum, functional diversity e functional redundancy
dens.convex<-(sum(matrix.resu[36,54],matrix.resu[36,59],matrix.resu[54,59]))/
  (factorial(3)/(factorial(2)*factorial(1)))
dens.convex
#estimadores e fisher alpha
dens.fisher<-(sum(matrix.resu[17,21]))/(factorial(2)/(factorial(2)*factorial(0)))
dens.fisher
mean.metric<-read.table("clipboard",header=TRUE)
rm(mean.metric)
