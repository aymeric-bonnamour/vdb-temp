#Ex 1
install.packages("ape")
library("ape")
??rtree
tree <- rtree(1000)
plot.phylo(tree)


#Ex 2
install.packages("leaflet")
library("leaflet")

m<-leaflet()
m<-addTiles(m, url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png")
m<-addMarkers(m, lng=4.85, lat=45.75, label="Ici Lyon")
m <- addMarkers(m, lng=2.35, lat=48.85, label="Paris")
m

#Load a new map from a data frame
loadMap <- function(dtframe){
  m <- leaflet(data=dtframe)
  m <- addTiles(m, url=" http://lifemap-ncbi.univ-lyon1.fr/osm_tiles/{z}/{x}/{y}.png", options=tileOptions(maxZoom = 42))
  return(m)
}
m <- load_map()
m

#Ex 3
ville <- c("Lyon", "Paris")
long <- c(4.85, 2.35)
lat <- c(45.75, 48.85)
table <- data.frame(ville, long, lat)
table
m <- leaflet(data=table)
m <- addTiles(m, url="https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png")
m <- addMarkers(m, lng=~long, lat=~lat, label=~ville)
m

#Ex 4
m<-leaflet()
m<-addTiles(m, url="http://lifemap-ncbi.univ-lyon1.fr/osm_tiles/{z}/{x}/{y}.png")
install.packages("jsonlite")
library(jsonlite)

GetCooFromTaxID<-function(taxids) {
  ##taxids is an array that contains taxids.
  ## url cannot be too long, so that we need to cut the taxids (100 max in one chunk)
  ## and make as many requests as necessary.
  taxids<-as.character(taxids) #change to characters.
  DATA<-NULL
  i<-1
  while(i<=length(taxids)) {
    cat(".")
    taxids_sub<-taxids[i:(i+99)]
    taxids_sub<-taxids_sub[!is.na(taxids_sub)]
    taxids_sub<-paste(taxids_sub, collapse="%20") #accepted space separator in url
    url<-paste("http://lifemap-ncbi.univ-lyon1.fr:8080/solr/taxo/select?q=taxid:(",taxids_sub,")&wt=json&rows=1000",sep="", collapse="")
    #do the request :
    data_sub<-fromJSON(url)
    DATA<-rbind(DATA,data_sub$response$docs[,c("taxid","lon","lat", "sci_name","zoom","nbdesc")])
    i<-i+100
  } 
  for (j in 1:ncol(DATA)) DATA[,j]<-unlist(DATA[,j])
  class(DATA$taxid)<-"character"
  return(DATA)
}

data<-GetCooFromTaxID(c(2,9443,2087))
data
m<-loadMap(data)
m <- addCircleMarkers(m, lng=~lon, lat=~lat)
m

##RÉCUPÉRER LES DONNÉES
EukGenomeInfo<-read.table("ftp://ftp.ncbi.nlm.nih.gov/genomes/GENOME_REPORTS/eukaryotes.txt", sep="\t", header=T, quote="\"", comment.char="")
## liste unique des taxid
taxids<-unique(EukGenomeInfo$TaxID)

## RÉCUPÉRER LES COORDONNÉES
DF<-GetCooFromTaxID(taxids)

## CALCULER LE NOMBRE DE GÉNOMES SÉQUENCÉS POUR CHAQUE TAXID
nbGenomeSequenced<-table(EukGenomeInfo$TaxID)
## l'ajouter à DF
DF$nbGenomeSequenced<-as.numeric(nbGenomeSequenced[DF$taxid])

##CALCULER LE NB DE GENOMES ENTIEREMENT ASSEMBLES POUR CHAQUE TAXID
##le calcul pour un seul taxid nommé 'tid' serait :
sum(EukGenomeInfo[which(EukGenomeInfo$TaxID=='tid'),]$Status=="Chromosome")
##on peut utiliser la fonction sapply pour le faire pour chaque taxid 
nbGenomeAssembled<-sapply(DF$taxid, function(x,tab) sum(tab[which(tab$TaxID==x),]$Status=="Chromosome"), tab=EukGenomeInfo)
DF$nbGenomeAssembled<-nbGenomeAssembled

##CALCULER LE TAUX de GC MOYEN  
tauxgcmoyen<-sapply(DF$taxid, function(x,tab) mean(as.numeric(as.character(tab[which(tab$TaxID==x),]$GC.)), na.rm=TRUE), tab=EukGenomeInfo)
DF$tauxgcmoyen<-tauxgcmoyen

##CALCULER LA TAILLE MOYENNE DES GÉNOMES EN Mb
SizeGenomeMb<-sapply(DF$taxid, function(x,tab) mean(tab[which(tab$TaxID==x),]$Size..Mb., na.rm=TRUE), tab=EukGenomeInfo)
DF$SizeGenomeMb<-SizeGenomeMb 
head(DF)

m <- loadMap(DF)
m <- addCircleMarkers(m, lng=~lon, lat=~lat, label=~sci_name)
m

##charger un package contenant des palettes de couleur
library(RColorBrewer) #ou "viridis" qui contient de bonnes couleurs aussi
## créer la fonction de palette
pal<-colorNumeric("Greens",0:100) ##green est une des palettes de RColorBrewer. Taper 
## display.brewer.all() pour les voir toutes

pal(12)