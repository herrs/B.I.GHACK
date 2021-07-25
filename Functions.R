address_root <<- "C:\\Users\\jenni\\Documents\\R\\Working Directory\\Hackathon"
api_key <- "5cfca50b1eb377dc77ab868e088438eeae09"
query <- "COVID*[Title/Abstract] AND english[LA] AND 2019[DP]"
limit = 0
paperlist <- paste(address_root, "\\basicsciencepapers.txt", sep="")
Mcol = PubMedSample$DE
max_words = 200
load_packages <- function(){
  a <- installed.packages()
  packages <- a[,1]
  pack_file = paste(address_root, "\\packagelist.txt",sep="")
  pack_list <- read.table(pack_file, sep=",")
  for(pack in pack_list)
  {if (!(is.element(pack, packages)))
  {install.packages(pack)}
   library(pack) 
  }
}
data_parse <- function(query, api_key, dbsource, format, limit){
  x <- pmQueryTotalCount(query = query, api_key = api_key)
  if(limit == 0)
  {limit = x$total_count}
  D <- pmApiRequest(query = query, limit = limit, api_key = api_key)
  M <- convert2df(file = D, dbsource = dbsource, format = format)
  M <- M %>% filter(DT == "JOURNAL ARTICLE"|DT =="LETTER")
  return(M)
}
PubMedData <- rbind(Past_DF, PubMedData)
# Modified biblioanalysis code from bibliometrix 
biblioAnalysis<-function(M,sep=";"){
  # initialize variables
  Authors=NULL
  Authors_frac=NULL
  FirstAuthors=NULL
  PY=NULL
  FAffiliation=NULL
  Affiliation=NULL
  Affiliation_frac=NULL
  CO=rep(NA,dim(M)[1])
  TC=NULL
  TCperYear=NULL
  SO=NULL
  Country=NULL
  DE=NULL
  ID=NULL
  MostCitedPapers=NULL
  
  
  
  
  
  # M is the bibliographic dataframe
  Tags<-names(M)
  
  if (!("SR" %in% Tags)){M=metaTagExtraction(M,"SR")}
  
  # temporal analyis
  
  #if ("PY" %in% Tags){Years=table(M$PY)}
  
  # Author's distribution
  
  if ("AF" %in% Tags){
    listAU=strsplit(as.character(M$AF),sep)
    listAU=lapply(listAU, function(l) trim(l))
    #nAU=unlist(lapply(listAU,length))  # num. of authors per paper
    nAU <- lengths(listAU)
    #fracAU=unlist(lapply(nAU,function(x){rep(1/x,x)}))  # fractional frequencies
    fracAU <- rep(1/nAU,nAU)
    AU=unlist(listAU)
    
    Authors=sort(table(AU),decreasing=TRUE)
    Authors_frac=aggregate(fracAU,by=list(AU),'sum')
    names(Authors_frac)=c("Author","Frequency")
    Authors_frac=Authors_frac[order(-Authors_frac$Frequency),]
    FirstAuthors=unlist(lapply(listAU,function(l){
      if (length(l)>0){l=l[[1]]} else {l=NA}
      return(l)
    }))
    
    AuSingleAuthoredArt=length(unique(FirstAuthors[nAU==1]))
    AuMultiAuthoredArt=length(Authors)-AuSingleAuthoredArt
  }
  
  #Total Citation Distribution
  if ("TC" %in% Tags){
    TC=as.numeric(M$TC)
    PY=as.numeric(M$PY)
    CurrentYear=as.numeric(format(Sys.Date(),"%Y"))
    TCperYear=TC/(CurrentYear-PY+1)
    if (!("DI" %in% names(M))) M$DI <- ""
    MostCitedPapers <- data.frame(M$SR,M$DI,TC,TCperYear,PY) %>%
      group_by(.data$PY) %>%
      mutate(NTC = .data$TC/mean(.data$TC)) %>%
      ungroup() %>% 
      select(-.data$PY) %>%
      arrange(desc(.data$TC)) %>%
      as.data.frame()
    
    names(MostCitedPapers)=c("Paper         ","DOI","TC","TCperYear","NTC")
  }
  
  # References
  nReferences <- 0
  if ("CR" %in% Tags){
    CR <- tableTag(M,"CR",sep)
    nReferences <- length(CR)
  }
  
  # ID Keywords
  if ("ID" %in% Tags){ID <- tableTag(M,"ID",sep)}
  
  # DE Keywords
  if ("DE" %in% Tags){DE=tableTag(M,"DE",sep)}
  
  # Sources
  if ("SO" %in% Tags){
    SO=gsub(",","",M$SO,fixed=TRUE)
    SO=sort(table(SO),decreasing = TRUE)
  }
  
  # All Affiliations, First Affiliation and Countries
  if (("C1" %in% Tags) & (sum(!is.na(M$C1))>0)){
    if(!("AU_UN" %in% Tags)){M=metaTagExtraction(M,Field="AU_UN")}
    AFF=M$AU_UN
    listAFF=strsplit(AFF,sep,fixed=TRUE)
    nAFF=unlist(lapply(listAFF,length))
    listAFF[nAFF==0]="NA"
    fracAFF=unlist(sapply(nAFF,function(x){
      if(x>0){x=rep(1/x,x)}else{
        x=0}
    }))  # fractional frequencies
    AFF=trim.leading(unlist(listAFF))  # delete spaces
    Affiliation=sort(table(AFF),decreasing=TRUE)
    Affiliation_frac=aggregate(fracAFF,by=list(AFF),'sum')
    names(Affiliation_frac)=c("Affiliation","Frequency")
    Affiliation_frac=Affiliation_frac[order(-Affiliation_frac$Frequency),]
    
    # First Affiliation
    FAffiliation=lapply(listAFF,function(l) l[1])
    
    
    
  }else{
    SCP_MCP=data.frame(Country=rep(NA,1),SCP=rep(NA,1))
  }
  if ("DT" %in% names(M)){
    Documents=table(M$DT)
    n=max(nchar(names(Documents)))
    names(Documents)=substr(paste(names(Documents),"                                              ",sep=""),1,n+5)
  }else{Documents=NA}
  
  results=list(Articles=dim(M)[1],             # Articles
               Authors=Authors,                # Authors' frequency distribution
               AuthorsFrac=Authors_frac,       # Authors' frequency distribution (fractionalized)
               FirstAuthors=FirstAuthors,      # First Author's list
               nAUperPaper=nAU,                # N. Authors per Paper
               Appearances=sum(nAU),            # Author appearances
               nAuthors=dim(Authors),          # N. of Authors
               AuMultiAuthoredArt=AuMultiAuthoredArt, # N. of Authors of multi-authored articles
               AuSingleAuthoredArt=AuSingleAuthoredArt, # N. of Authors of single-authored articles
               MostCitedPapers=MostCitedPapers,# Papers sorted by citations
               Years=PY,                       # Years
               FirstAffiliation=unlist(FAffiliation),  # Affiliation of First Author
               Affiliations=Affiliation,       # Affiliations of all authors
               Aff_frac=Affiliation_frac,      # Affiliations of all authors (fractionalized)
               TotalCitation=TC,               # Total Citations
               TCperYear=TCperYear,            # Total Citations per year
               Sources=SO,                     # Sources
               DE=DE,                          # Keywords
               ID=ID,                          # Authors' keywords
               Documents=Documents,
               nReferences = nReferences,      # N. of References
               DB=M$DB[1])
  class(results)<-"bibliometrix"
  
  return(results)
}
topk_paper_summary <- function(M, k)
{
results <- biblioAnalysis(M)
x <- summary(results)
capture.output(x, file= paste(address_root, "basicsciencesummary.txt", sep=""))
write(results$MostCitedPapers$DOI[1:k], file= paste(address_root, "basicsciencepapers.txt", sep=""), ncolumns = k, sep="|")
write(x$MostProdAuthors$`Authors       `, file= paste(address_root, "mostprodauthors.txt", sep=""),ncolumns = k, sep="|")
}

get_landmark_papers <- function(M, paperlist)
{LandmarkPapers <- read.table(paperlist)
LM <- M %>% filter(grepl(LandmarkPapers, DI))
return(LM)}

plot_country <- function(M)
{CountryOfOrigin <- ggplot(data=M, aes(x=factor(SO_CO), fill=factor(SO_CO)))+ geom_bar(stat="count", width=0.7) + scale_fill_brewer(palette="Dark2")+ 
  scale_colour_discrete(name = "Whatever I please")
return(CountryOfOrigin)
}

word_cloud <- function(Mcol, max_words)
{
write.table(Mcol, file = paste(address_root, "\\wordsoup.txt", sep=""), sep = ",", eol = ", ", row.names = FALSE)
text <- readLines(paste(address_root, "\\wordsoup.txt"))
text <- strsplit(text, ";")
docs <- Corpus(VectorSource(text))
toSpace <- content_transformer(function(x, pattern) gsub(pattern," ",x))

docs <- tm_map(docs, content_transformer(tolower))
# Remove numbers
docs <- tm_map(docs, removeNumbers)
# Remove english common stopwords
docs <- tm_map(docs, removeWords, stopwords("english"))
# Remove punctuations
docs <- tm_map(docs, removePunctuation)
# Remove obvious keywords
docs <- tm_map(docs, removeWords, c("health", "covid", "coronavirus", "pandemic", "sarscov", "disease", "research", "clinical"))
# Eliminate extra white spaces
docs <- tm_map(docs, stripWhitespace)
dtm <- TermDocumentMatrix(docs)
m <- as.matrix(dtm)
v <- sort(rowSums(m),decreasing=TRUE)
d = data.frame(word = names(v),freq=v)

set.seed(1234)
cloud <- wordcloud(words = d$word, freq = d$freq, min.freq = 1, max.words=max_words, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
return(cloud)
}
network_build <- function(M)
{TopIAuthors <- readLines(paste(address_root, "mostprodauthors.txt", sep=""))
PubMedTopAuthors <- M %>% filter(grepl(TopIAuthors, AF))


#run modified biblioNetwork script 
NetMatrix <-biblioNetwork(PubMedTopAuthors, analysis = "collaboration", network = "authors", sep=";", type="sparse")
#NetMatrix <- cocMatrix(PubMedData, Field = "AF", type="sparse", sep=";")
net=networkPlot(NetMatrix, n=150, Title = "Institution collaboration", type="auto",size = 10, size.cex=T, edgesize=3,labelsize=0.01)

#######
#NODE DEGREES AND CENTRALITY
stats <- networkStat(NetMatrix, stat="all") #use for centrality
NodeCentralities <- data.frame(stats$vertex$vertexID, stats$vertex$vertexCentrDegree)
attach(NodeCentralities)
NodeCentralities <- NodeCentralities[order(-stats.vertex.vertexCentrDegree), ]
detach(NodeCentralities)

NodeDegree <- degree(NetMatrix)
NodeDegree <- data.frame(stats$vertex$vertexID, NodeDegree)#use for node degrees

NetMatrix <- as.matrix(NetMatrix)
PubMedTopAuthors.statnet = as.network(NetMatrix, stat="all",directed = FALSE, names.eval = "edge.lwd", ignore.eval = FALSE)
PubMedTopAuthors.statnet # view network summary

plot.network(PubMedTopAuthors.statnet, edge.col = "gray", edge.lwd = "edge.lwd", label = "vertex.names", label.cex = .5, label.pad = 0, label.pos = 1)

PubMedTopAuthors.nodes <- data.frame(id = 1:length(PubMedTopAuthors.statnet%v%"vertex.names"),
                                     label = PubMedTopAuthors.statnet%v%"vertex.names",
                                     title = PubMedTopAuthors.statnet%v%"vertex.names",
                                     size = 5*(2+PubMedTopAuthors.statnet%v%"size"))

PubMedTopAuthors.edges <- data.frame(from=data.frame(as.edgelist(PubMedTopAuthors.statnet))$X1, 
                                     to=data.frame(as.edgelist(PubMedTopAuthors.statnet))$X2)
PubMedTopAuthors_interactive = visNetwork(PubMedTopAuthors.nodes, PubMedTopAuthors.edges, main = "Top Author Co-Author Network", width = 800, height = 800) %>% 
  visIgraphLayout(layout = "layout_nicely", type = "full")
PubMedTopAuthors_interactive = PubMedTopAuthors_interactive  %>%
  visNodes(color = list(background = "white", highlight = "red", hover = list(border = "red"))) %>%
  visEdges(selectionWidth = 10, color = list(highlight = "#2B7CE9")) # view interactive network
PubMedTopAuthors_interactive = PubMedTopAuthors_interactive  %>%  
  visOptions(nodesIdSelection = list(enabled  = TRUE, useLabels = TRUE, main = "Select by Author"))
return(PubMedTopAuthors_interactive)
}
freqToVec <- function(df){
  x <- replicate(df[2], df[1])
  return(x)
  
}

reduceWordout <- function(words, df, i, full){
  newDF <- df
  reduce <- reduceWordin(words[i], newDF)
  
  
  full <- data.frame(Var1 = words[i] , Freq = sum(reduce$Freq))
  newDF <- reduceDF(newDF, reduce,i)
  newDF <- rbind(newDF, full)
  i = i +1
  len = length(words)
  if (i <= len){
    reduceWordout(words, newDF, i, full)
    
  }else{
    return(newDF)
  } 
}

reduceWordin <- function(word, df){
  reduce <- df[grep(word, df$Var1),]
  
  return(reduce)  
}

reduceDF <- function(reduced, subset,i){
  df <- anti_join(reduced, subset)
  
  
  return(df)
  
}

getDepartments <- function(fullDF){
  
  ##filter csv file t get department column
  #take the first string before comma
  
  fullDF <- fullDF %>% filter(DT == "JOURNAL ARTICLE" )
  fullDF <- fullDF %>% select( C1)
  fullDF <- sapply(fullDF, strsplit, ";")
  fullDF <- unlist(sapply(fullDF,data.frame ))
  fullDF <- data.frame(substr(fullDF,1,regexpr(",",fullDF)-1))
  fullDF <- sapply(fullDF, tolower)
  fullDF <- as.vector((fullDF))
  
  
  
  
  ##take words after for, and, of
  department <- sub(".*of" ,"", fullDF)
  department <- sub(".*and " ,"", department)
  department <- sub(".*for " ,"", department)
  
  department <- as.data.frame(table(unlist(department)))
  department <- with(department, department[ Var1 != "", ] )
  department$Var1 <- as.character(department$Var1)
  department <- department[order(department$Freq, decreasing = T), ]
  
  
  
  
  
  ###biblioAnalysis stuff
  
  
  
  
  
  ################################################# ATTEMPTS to collapse rows with medicine
  #in it to a single row with updated counts
  
  
  
  full <- data.frame(field = as.character(), frequency = as.numeric())
  
  words <- c( "psy\\w*", "stat\\w*", "viro\\w+", "child",  "card\\w+" , "imaging",  "neuro\\w+", "radio\\w+", "stroke", "bioinformatic\\w*", "behav\\w*", "econ\\w*",  "mech\\w*",  "disease\\w*", "data\\w*",  "epidemiology","policy\\w*", "pedia\\w*|paed\\w*",  "immun\\w*",  "comp\\w*", "pharm\\w",  "physics\\w*", "math\\w*", "agri\\w*", "geno\\w*|gene\\w*",  "public health", "immunology",  "brain",  "biolog\\w*","chem\\w*", "engi\\w*",  "med*|hosp\\w*",  "surg\\w*", "model\\w*" )
  
  
  
  reduce <- reduceWordout(words, department,1, full)
  reduce <- reduce[order(reduce$Freq, decreasing = T),]
  reduce <- reduce[reduce$Freq >3,]
  nice <- reduce[1:20,] 
  
  ##name it nicely
  nice$Var1 <- c("Medicine", " Neurology", "Disease", "Biology", "Cardiology", "Surgery", "Psychology", "Radiology", "Pediatrics", "Pharmaceutical", "Public Health", "Genetics", "Immunology", "Stroke Research", "Virology", "Statistics", "Engineering", "Chemistry", "Computer Science", "Epidemology")
  names(nice) <- nice[1,]
  new <- apply(data.frame(nice), 1, freqToVec)
  new <- data.frame(unlist(new))
  names(new) <- "field"
  
  
  
  pic <-  ggplot(data.frame(new),aes(x=reorder(field,field,
                                               function(x)-length(x)))) + geom_bar(aes(fill = field)) +xlab("Field") + theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.2)) + theme(legend.position = 'none')
  
}
