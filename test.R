#Z-test for difference between mutational and evolutionary correlations
#Use setwd() first

library(phytools)
library(geiger)
library(mvtnorm)

d1<-read.table("list_strain.txt",header=TRUE,sep="\t");l1=as.character(d1[,1]);clade=as.character(d1[,2])
d2<-read.table("isolate_list_new.txt",header=TRUE,sep="\t");l2=as.character(d2[,1]);l3=as.character(d2[,2])
dm<-read.table("1011DistanceMatrixBasedOnSNPs.tab",header=TRUE,sep="\t")
rownames(dm)=as.character(dm[,1]);dm=dm[,2:ncol(dm)]
dm=data.matrix(dm)
dm=(dm/100)*1544489
#dm=(dm/100)*1.21e7
#dm=dm*1000

#generate a list of strains that belong to unambiguous clades and are available in all datasets
rm=c()
for(i in 1:length(l1)){
	check=which((l2)==l1[i])
	if(length(check)==0|(clade[i]=="none")){
		rm=c(rm,i)
	}
}
l1=l1[-rm] #retain overlap

list=c()
for(i in 1:nrow(dm)){
	name=l2[which((l3)==rownames(dm)[i])]
	if(length(which((l1)==name))>0){
		list=c(list,i)
		rownames(dm)[i]=name
		colnames(dm)[i]=name
	}
}
dm=dm[list,list]

tr<-nj(dm)
tr=multi2di(tr)
obs=c()
for(i in 1:nrow(tr$edge)){
	descend=which(tr$edge[,1]==tr$edge[i,2])
	if(length(descend)==0){
		obs=c(obs,i)
	}
}
sp=tr$edge[obs,]

wt<-read.table("wt_mean.txt",sep="\t")
wt=colMeans(wt)
trait_del<-read.table("220_list.txt",sep="\t")
module=as.character(trait_del[,3])
trait_del=as.character(trait_del[,2])
del<-read.table("mt_mean.txt",sep="\t")
del_norm=matrix(0,nrow=nrow(del),ncol=ncol(del))
for(i in 1:ncol(del)){
	del_norm[,i]=log(del[,i]/wt[i])
}

colnames(del_norm)=trait_del
ntrait=ncol(del_norm);npair=(ntrait^2-ntrait)/2
mmat=cov(del_norm);mcor=cov2cor(mmat)
em=eigen(mcor)$values

div<-read.table("div_mean_combined.txt",sep="\t")
div<-div[1:37,]
div=div[-rm,]
div_norm=matrix(0,nrow=nrow(div),ncol=ncol(div))
for(i in 1:ncol(div)){
	div_norm[,i]=log(div[,i]/wt[i])
}
colnames(div_norm)=trait_del
rownames(div_norm)=l1
rmat=ratematrix(tr,div_norm)
rcor=cov2cor(rmat)
er=eigen(rcor)$values

pcor_m=matrix(0,nrow=ntrait,ncol=ntrait)
pcor_r=matrix(0,nrow=ntrait,ncol=ntrait)
for(i in 1:ntrait){
	for(j in 1:ntrait){
		if(i<j){
			tvm=abs(sqrt(nrow(del_norm)-2)*mcor[i,j]/sqrt(1-(mcor[i,j]^2)))
			tvr=abs(sqrt(nrow(sp)-2)*rcor[i,j]/sqrt(1-rcor[i,j]^2))
			pcor_m[i,j]=2*pt(tvm,df=nrow(del_norm)-2,lower.tail=FALSE)
			pcor_r[i,j]=2*pt(tvr,df=nrow(sp)-1-2,lower.tail=FALSE)
			pcor_m[j,i]=pcor_m[i,j]
			pcor_r[j,i]=pcor_r[i,j]
		}
	}
}

#direct comparison of two correlation matrices
trait1=c();trait2=c();cor1=c();cor2=c();pcor1=c();pcor2=c()
pv=c()
for(i in 1:ntrait){
	for(j in 1:ntrait){
		if(i<j){
			trait1=c(trait1,colnames(div_norm)[i])
			trait2=c(trait2,colnames(div_norm)[j])
			cor1=c(cor1,mcor[i,j])
			cor2=c(cor2,rcor[i,j])
			pcor1=c(pcor1,pcor_m[i,j])
			pcor2=c(pcor2,pcor_r[i,j])
			z1=0.5*(log(1+mcor[i,j])-log(1-mcor[i,j]))
			z2=0.5*(log(1+rcor[i,j])-log(1-rcor[i,j]))
			z=(z2-z1)/sqrt((nrow(del_norm)-3)^-1+(nrow(sp)-1-3)^-1)
			pv=c(pv,2*pnorm(-abs(z)))
		}
	}
}
padj=p.adjust(pv,method="fdr")
cor_inf=data.frame(trait1,trait2,cor1,pcor1,cor2,pcor2,cor2/cor1,pv,padj)
colnames(cor_inf)=c("trait1","trait2","COR_M","P_M","COR_E","P_E","ratio","P","P_adjust")

#simulation following BM
Ntest=1000
cor_all=list()
pv_all=list()
padj_all=list()
type_all=list()
nsignif=c()
nsignif_adj=c()
nup=c();ndn=c();nrv=c()
for(test in 1:Ntest){
	#simulation
	phe=matrix(0,nrow=nrow(tr$edge),ncol=ntrait)
	for(i in 1:nrow(tr$edge)){
		ances=which(tr$edge[,2]==tr$edge[i,1])
		if(length(ances)==0){
			phe_ances=rep(0,ntrait)
		}else{
			phe_ances=phe[ances,]
		}
		phe[i,]=phe_ances+tr$edge.length[i]*rmvnorm(1,sigma=mmat)
	}
	phe_obs=phe[obs,]
	rownames(phe_obs)=rep(0,nrow(phe_obs))
	for(i in 1:nrow(phe_obs)){
		rownames(phe_obs)[i]=tr$tip[sp[i,2]]
	}
	rmat_control=ratematrix(tr,phe_obs)
	rcor_control=cov2cor(rmat_control)
	cor_all[[test]]=rcor_control

	#test for pairwise correlations
	pv_vect=c() #store P-values in a vector
	type_vect=c()
	for(i in 1:ntrait){
		for(j in 1:ntrait){
			if(i<j){
				z1=0.5*(log(1+mcor[i,j])-log(1-mcor[i,j]))
				z2=0.5*(log(1+rcor_control[i,j])-log(1-rcor_control[i,j]))
				z=(z2-z1)/sqrt((nrow(del_norm)-3)^-1+(nrow(sp)-1-3)^-1)
				pv_vect=c(pv_vect,2*pnorm(-abs(z)))
				#type of difference
				type_pair="none"
				if(rcor_control[i,j]*mcor[i,j]<0){
					if((pcor_m[i,j]<0.05)&(pcor_r[i,j]<0.05)){
						type_pair="rv"
					}else{
						if(pcor_m[i,j]>0.05){
							type_pair="up"
						}
						if(pcor_r[i,j]>0.05){
							type_pair="dn"
						}
					}
				}else{
					if(abs(rcor_control[i,j])>abs(mcor[i,j])){
						type_pair="up"
					}else{
						type_pair="dn"
					}
				}
				type_vect=c(type_vect,type_pair)
			}
		}
	}
	pv_all[[test]]=pv_vect
	nsignif=c(nsignif,length(which(pv_all[[test]]<0.05)))
	padj_all[[test]]=p.adjust(pv_vect,method="fdr")
	nsignif_adj=c(nsignif_adj,length(which(padj_all[[test]]<0.05)))
	type_vect[which(padj_all[[test]]>0.05)]="none"
	type_all[[test]]=type_vect
	nup=c(nup,length(which(type_all[[test]]=="up")))
	ndn=c(ndn,length(which(type_all[[test]]=="dn")))
	nrv=c(nrv,length(which(type_all[[test]]=="rv")))
}
median(nsignif_adj)
median(nup);median(ndn);median(nrv)


#for each trait pair, check in how many simulations there's significant difference
num=c();num_adj=c()
for(i in 1:nrow(cor_inf)){
	pv_pw=c();padj_pw=c()
	for(n in 1:Ntest){
		pv_pw=c(pv_pw,pv_all[[n]][i])
		padj_pw=c(padj_pw,padj_all[[n]][i])
	}
	num=c(num,length(which(pv_pw<0.05)))
	num_adj=c(num_adj,length(which(padj_pw<0.05)))
}

#mean correlation across simulations for each trait pair
sim_mean=c()
for(i in 1:ntrait){
	for(j in 1:ntrait){
		if(i<j){
			cor.est=c()
			for(n in 1:Ntest){
				cor.est=c(cor.est,cor_all[[n]][i,j])
			}
			sim_mean=c(sim_mean,mean(cor.est))
		}
	}
}

out=data.frame(cor_inf,num,num_adj,sim_mean)
write.table(out,file="cor_inf_yeast.txt",sep="\t")
write.table(nsignif,file="nsignif.txt",sep="\t")
write.table(data.frame(nsignif_adj,nup,ndn,nrv),file="nsignif_adj_yeast.txt",sep="\t")

length(which(cor_inf$ratio>1&cor_inf$P_adjust<0.05))+length(which(cor_inf$ratio<0&cor_inf$P_M>0.05&cor_inf$P_adjust<0.05))
length(which(cor_inf$ratio>0&cor_inf$ratio<1&cor_inf$P_adjust<0.05))+length(which(cor_inf$ratio<0&cor_inf$P_E>0.05&cor_inf$P_adjust<0.05))
length(which(cor_inf$ratio<0&cor_inf$P_adjust<0.05&cor_inf$P_M<0.05&cor_inf$P_E<0.05))

#integration/modularity analysis
vobs=var(er)
Ntest=5000
vnull=rep(0,Ntest);rknull=rep(0,Ntest);pfk=rep(0,Ntest)
for(test in 1:Ntest){
	subset=del_norm[sample(1:ntrait,nrow(sp)-1),]
	enull=eigen(cov2cor(cov(subset)))$values
	vnull[test]=var(enull);rknull=qr(cov2cor(cov(subset)))$rank
	pfk[test]=fligner.test(list(enull,er))$p.value
}
tail=length(which(vnull>vobs))
2*min(c(tail,length(vnull)-tail))/length(vnull) #2-side p-value

cr_cal<-function(m,mod){
	m11=m[which((module)==mod),which((module)==mod)]
	m22=m[which((module)!=mod),which((module)!=mod)]
	m12=m[which((module)==mod),which((module)!=mod)]
	m21=m[which((module)!=mod),which((module)==mod)]
	m11a=m11;diag(m11a)=0
	m22a=m22;diag(m22a)=0
	numerator=sum(diag(m12%*%m21))
	denomenator=sqrt(sum(diag(m11a%*%m11a))*sum(diag(m22a%*%m22a)))
	cr=sqrt(numerator/denomenator)
	return(cr)
}

cr_obs=mean(c(cr_cal(rmat,"actin"),cr_cal(rmat,"cell wall"),cr_cal(rmat,"nucleus")))
cr_null=rep(0,Ntest)
for(test in 1:Ntest){
	subset=del_norm[sample(1:ntrait,nrow(sp)),]
	rmat_control=cov(subset)
	cr_null[test]=mean(c(cr_cal(rmat_control,"actin"),cr_cal(rmat_control,"cell wall"),cr_cal(rmat_control,"nucleus")))
}
tail=length(which(cr_null>cr_obs))

