### script - data analysis adpated from Ribas Deulofeu et al. ###
# update May 9th

rm(list=ls()) # cleaning memory 

## Working directory (to customize)
# setwd('C:/Users/Vdenis/OneDrive/Documents/Class/OCEAN_7194_FIELD/Lectures/Part 6') # set up working directory

## Packages and functions
library(vegan) # install.packages('vegan')
library(gclus) # install.packages('gclus')
source('Function/coldiss.R')


## Dataset and factors
### Import dataset and remove unstable substrate (transect level)
data<- read.csv('Data/TWdataset.csv',h=TRUE,row.names='OTUs') # import dataset
str(data)
data1<-as.data.frame(data) # conversion as a dataframe
data1[,'sand']<- NULL    #  delete unstable substrate 'sand'
data1[,'gravel_small_rubble']<- NULL #  delete unstable substrate 'gravel_small_rubble'
data1[,'rubble']<- NULL #  delete unstable substrate 'gravel_small_rubble'
data1<-as.matrix(data1) # conversion as a matrix
data2<-prop.table(data1,1) # conversion in percentage
### Import factor
fact1<-read.csv('Data/TWfactors.csv', h=TRUE,row.names=1)
fact1

## Convert dataset to site level
SiteData<- aggregate(data2, by=list(Category=fact1$site), FUN=mean) # mean OTUs percentages per sites
rownames(SiteData)<-SiteData$Category # apply Category (Site) as rownames
SiteData[,'Category']<- NULL    #  delete 'Category' column as variable 

# factor: region (site level) #
SiteName<- c('Bitou','Chaikou','Chinwan_Inner_Bay', 'Cimei', 'Dabaisha','Dingbaisha',  'Gongguan', 'Gupoyu','Hongchai', 'Houbihu', 'Jialeshuei', 'Keelung_Island','Leidashih', 'Longdon','Longkeng', 'Outlet', 'Pon_Pon_Tan','Sangjiaowan', 'Shihland', 'Siyuping','Tiaoshih',  'Tanzihwan','Wa_En_Tung', 'Wanlitung','Yeliu')
RegionName<- c('North_Taiwan','Green_Island','Penghu','Penghu','Green_Island','Kenting','Green_Island','Penghu','Kenting','Kenting','Kenting','North_Taiwan','Kenting','North_Taiwan','Kenting','Kenting','Penghu','Kenting','Green_Island','Penghu','Kenting','Kenting','Penghu','Kenting','North_Taiwan')
fact1<-data.frame(SiteName,RegionName) # regional factor at site level


## distance & cluster##
dis1<- vegdist(SiteData, method='bray') # bray-curtis similarity (site level)
coldiss(dis1,byrank=F,diag=T) # for the bc dissimilarity on raw data 
clus1<-hclust(dis1, method='average') # cluster average linkage
plot(clus1)

## nonmetric Multidimensional Scaling ##         
mds1site <- metaMDS(SiteData, dist='bray', try=99,type='n') # nonmetric Multidimensional Scaling analysis 


plot(mds1site, display='sites', type='none') # plot Nonmetric Multidimensional Scaling



text(mds1site, 'species',col = rgb(0,0,0,alpha=0.3) , cex=0.5) # visualize species names



points(mds1site, 'sites', pch=19,cex=2, col='yellow',select=fact1$RegionName=='Penghu') # add transects from Penghu
points(mds1site, 'sites', pch=19, cex=2, col='green',select=fact1$RegionName=='Green_Island')# add transects from Green Island
points(mds1site, 'sites', pch=19, cex=2,col='cyan',select=fact1$RegionName=='North_Taiwan') # add transects from North Taiwan
points(mds1site, 'sites', pch=19, cex=2,col='red',select=fact1$RegionName=='Kenting') # add transects from Kenting
ordispider (mds1site, RegionName, col='red') # position centroids

## similarity percentages ##
simper1site<-simper(SiteData, RegionName) # similarity percentages analysis
simper1site # discriminating species between regions

## dispersion test ## 
mod<-betadisper (dis1,RegionName) # multivariate homogeneity of groups dispersions on untransformed data (site Level)  
permutest (mod, pairwise=T) # ANOVA like permutation test (see ?permutest)
TukeyHSD(mod) # pairwise comparison
boxplot(mod) # visualization of the dispersion


## PERMANOVA ##
permanovaRegion<-adonis2(formula=SiteData~RegionName, data=fact1, permutations=9999, method='bray') # Permutational Multivariate Analysis of Variance using bray-curtis  
permanovaRegion # permanova results

## ANOSIM ##

anosim1<-anosim(dis1,RegionName, permutations=9999) # analysis of similarities
anosim1 # results anosim

