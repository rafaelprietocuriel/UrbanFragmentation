
#### code designed to construct metrics of fragmentation
#### for keeping track of which cities have been computed, it creates a folder and files
#### in /HasRun it will create an empty DF with the name of cities that have been computed
#### each time it starts, it runs through /HasRun to see what has been computed

#### inputs
#### uses Africapolis in a csv format
#### and the buildings dataset, already in the following form:
#### for each city, we identify the buildings inside the polygon
#### then only the coordinates and area of the buildings of that city will be used

#### outputs
#### it will create many folders to manage the separate type of files
#### the outputs include:
#### a set of metrics for each urban form
#### that are added to the Africapolis database as variables
#### a polycentrism tree
#### a density figure in pdf and RData
#### an analysis of the size of buildings

#### packages
{
    require(geosphere)
    require(pracma)
    require(spatstat)
    require(maptools)
    require(denpro)
    memory.limit(size = 400000)
    options(try.outFile = stdout()) 
}

#### create outputs dirs
{
    dir.create("polTrees")  
    dir.create("DensityImages")
    dir.create("AfOutput")
    dir.create("figures")
    dir.create("BuildingSize")
    dir.create("HasRun")
}

#### functions
{
    PreClean <- function(S){
        if(!is.numeric(S$area_in_meters)){
            S <- S[-which(S$area_in_meters == "area_in_meters"), ]
            S$area_in_meters <- as.numeric(levels(S$area_in_meters))[S$area_in_meters]
            S$confidence <- as.numeric(levels(S$confidence))[S$confidence]
            S$lon <- as.numeric(levels(S$lon))[S$lon]
            S$lat <- as.numeric(levels(S$lat))[S$lat]
        }
        names(S) <- c("area", "confidence", "agglosID", "x", "y")
        return(S)
    }
    
    FragmentationI <- function(S){
        ST <- S #backup
        if(dim(S)[1] > 10000){
            u <- sample.int(dim(S)[1], size = 10000)
            S <- S[u, ]}
        D <- distm(data.frame(S$x, S$y))/1000
        diag(D) <- NA
        P <- matrix(rep(S$area, dim(S)[1]), ncol = dim(S)[1])
        W <- P * t(P)
        diag(W) <- NA
        W <- W/sum(W, na.rm = T)
        ### for yellow + green
        S <- ST[ST$confidence > 0.65, ]
        if(dim(S)[1] > 10000){
            u <- sample.int(dim(S)[1], size = 10000)
            S <- S[u, ]}
        DY <- distm(data.frame(S$x, S$y))/1000
        diag(DY) <- NA
        P <- matrix(rep(S$area, dim(S)[1]), ncol = dim(S)[1])
        WY <- P * t(P)
        diag(WY) <- NA
        WY <- WY/sum(WY, na.rm = T)
        ### for only green
        S <- ST[ST$confidence > 0.7, ]
        if(dim(S)[1] > 10000){
            u <- sample.int(dim(S)[1], size = 10000)
            S <- S[u, ]}
        DG <- distm(data.frame(S$x, S$y))/1000
        diag(DG) <- NA
        P <- matrix(rep(S$area, dim(S)[1]), ncol = dim(S)[1])
        WG <- P * t(P)
        diag(WG) <- NA
        WG <- WG/sum(WG, na.rm = T)
        
        return(list(MeanDistanceBuildings = mean(D, na.rm = T),
                    WMeanDistanceBuildings = sum(W*D, na.rm = T),
                    MaxBuildingsDistance = max(D, na.rm = T),
                    MeanDistanceBuildingsY = mean(DY, na.rm = T),
                    WMeanDistanceBuildingsY = sum(WY*DY, na.rm = T),
                    MaxBuildingsDistanceY = max(DY, na.rm = T),
                    MeanDistanceBuildingsG = mean(DG, na.rm = T),
                    WMeanDistanceBuildingsG = sum(WG*DG, na.rm = T),
                    MaxBuildingsDistanceG = max(DG, na.rm = T)))
    }
    
    PolycentricI <- function(S, bw = 0.001){ #bw = 0.005 better
        pcf<-pcf.kern(as.matrix(cbind(S$x, S$y)),
                      weights = S$area,
                      h= bw,
                      N=c(200,200))
        lst<-leafsfirst(pcf)
        td<-treedisc(lst,pcf,ngrid=100)
        levs <- td$level
        vol <- td$volume
        netTree <- td$parent
        
        #### compute base
        Base <- 0
        Parts <- c()
        u <- which.max(levs)
        NodePar <- netTree[u]
        Base <- Base + vol[u]
        levs[u] <- 0
        while(NodePar > 0){
            u <- NodePar
            NodePar <- netTree[u]
            Base <- Base + vol[u]
            levs[u] <- 0
            vol[u] <- 0
        }
        
        #### compute branches
        while(sum(levs)>0){
            u <- which.max(levs)
            NodePar <- netTree[u]
            tempBase <- vol[u]
            levs[u] <- 0
            vol[u] <- 0
            while(NodePar > 0){
                u <- NodePar
                NodePar <- netTree[u]
                tempBase <- tempBase + vol[u]
                levs[u] <- 0
                vol[u] <- 0
            }
            Parts <- c(Parts, tempBase)
            
        }
        WeighBranch <- c(Base, Parts)/Base
        return(list(NBranches = length(WeighBranch),
                    Polycentrism = sum(WeighBranch*(1:length(WeighBranch))),
                    HerfindahlIndex = sum((WeighBranch/sum(WeighBranch))^2),
                    netT = td))
    }
    
    DensityI <- function(S, bw = 0.001){
        CityEdge <- owin(poly=list(x=c(min(S$x), max(S$x), max(S$x), min(S$x)),
                                   y=c(min(S$y), min(S$y), max(S$y), max(S$y))))
        Buildings <- ppp(S$x, S$y,
                         marks = S$area,
                         window = CityEdge)
        D <- density.ppp(Buildings, 
                         sigma = bw,
                         weights = Buildings$marks)
        mx <- which(D$v == max(D$v), arr.ind = T)
        CentreLon <- D$xcol[mx[2]]
        CentreLat <- D$yrow[mx[1]]
        CentreDensity <- max(D$v)
        DC <- distm(data.frame(S$x, S$y),
                    data.frame(CentreLon, CentreLat))
        NBuildingsCentre1km <- sum(DC < 1000)
        TotalFootprintCentre1km <- sum(S$area[DC < 1000])
        NBuildingsCentre3km <- sum(DC < 3000)
        TotalFootprintCentre3km <- sum(S$area[DC < 3000])
        MeanArea0100m <- mean(S$area[DC < 0100 ])
        MeanArea0250m <- mean(S$area[DC < 0250 ])
        MeanArea0500m <- mean(S$area[DC < 0500 ])
        MeanArea0750m <- mean(S$area[DC < 0750 ])
        MeanArea1000m <- mean(S$area[DC < 1000 ])
        MeanArea1500m <- mean(S$area[DC < 1500 ])
        MeanArea2000m <- mean(S$area[DC < 2000 ])
        MeanArea2500m <- mean(S$area[DC < 2500 ])
        MeanArea3000m <- mean(S$area[DC < 3000 ])
        MeanArea4000m <- mean(S$area[DC < 4000 ])
        
        return(list(CentreLon =CentreLon, 
                    CentreLat=CentreLat, 
                    CentreDensity=CentreDensity, 
                    NBuildingsCentre1km=NBuildingsCentre1km,
                    TotalFootprintCentre1km=TotalFootprintCentre1km,
                    NBuildingsCentre3km=NBuildingsCentre3km,
                    TotalFootprintCentre3km=TotalFootprintCentre3km,
                    MeanArea0100m = MeanArea0100m,
                    MeanArea0250m = MeanArea0250m,
                    MeanArea0500m = MeanArea0500m,
                    MeanArea0750m = MeanArea0750m,
                    MeanArea1000m = MeanArea1000m,
                    MeanArea1500m = MeanArea1500m,
                    MeanArea2000m = MeanArea2000m,
                    MeanArea2500m = MeanArea2500m,
                    MeanArea3000m = MeanArea3000m,
                    MeanArea4000m = MeanArea4000m,
                    DensityImage = D))
            
    }
    
    AreaDistDecreasing <- function(S){
        x <-  rev(c(0, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200, 250, 300, 400, 500, 750, 1000, 1500, 2500, 5000, 7500, 10000, 20000, 50000, 100000,1000000))  
        nBuildings <- 0*x
        aBuildings <- 0*x
        for (k in 1:length(x)){
            ind <- S$area > x[k]
            nBuildings[k] <- sum(ind)
            aBuildings[k] <- sum(S$area[ind])
        }
            return(data.frame(buildingSize = x,
                          nBuildings = nBuildings,# / dim(S)[1],
                          aBuildings = aBuildings))# / sum(S$area_in_meters)))
    }
}

#### input data Africapolis
#### uses Africapolis to reac city by city polygons
{
Af <- read.csv("Africapolis.csv", encoding="UTF-8")
names(Af)[1] <- "agglosID"
names(Af)[4] <- "x"
names(Af)[5] <- "y"
n <- dim(Af)[1]
Af <- Af[order(Af$agglosID, decreasing = F),]
#### create all variables and set them with 0 value
Af$NBuildings <- rep(0, n)
Af$TotalFootprint <- rep(0, n)
Af$AreaConfidence <- rep(0, n)
Af$MeanBuildingArea <- rep(0, n)
Af$MedianBuildingArea <- rep(0, n)
Af$VarBuildingArea <- rep(0, n)
Af$SmallBuildings <- rep(0, n)
Af$LargeBuildings <- rep(0, n)
Af$AreaBuilt <- rep(0, n)
Af$NBranches <- rep(0, n)
Af$Polycentrism <- rep(0, n)
Af$HerfindahlIndex <- rep(0, n)
Af$MeanDistanceBuildings <- rep(0, n)
Af$WMeanDistanceBuildings <- rep(0, n)
Af$MaxBuildingsDistance <- rep(0, n)
Af$MeanDistanceBuildingsY <- rep(0, n)
Af$WMeanDistanceBuildingsY <- rep(0, n)
Af$MaxBuildingsDistanceY <- rep(0, n)
Af$MeanDistanceBuildingsG <- rep(0, n)
Af$WMeanDistanceBuildingsG <- rep(0, n)
Af$MaxBuildingsDistanceG <- rep(0, n)
Af$CentreLon <- rep(0, n)
Af$CentreLat <- rep(0, n)
Af$CentreDensity <- rep(0, n)
Af$NBuildingsCentre1km <- rep(0, n)
Af$TotalFootprintCentre1km <- rep(0, n)
Af$NBuildingsCentre3km <- rep(0, n)
Af$TotalFootprintCentre3km <- rep(0, n)
Af$MeanArea0100m <- rep(0, n)
Af$MeanArea0250m <- rep(0, n)
Af$MeanArea0500m <- rep(0, n)
Af$MeanArea0750m <- rep(0, n)
Af$MeanArea1000m <- rep(0, n)
Af$MeanArea1500m <- rep(0, n)
Af$MeanArea2000m <- rep(0, n)
Af$MeanArea2500m <- rep(0, n)
Af$MeanArea3000m <- rep(0, n)
Af$MeanArea4000m <- rep(0, n)
Af$FootpLargeBuildings050m <- rep(0, n)
Af$FootpLargeBuildings100m <- rep(0, n)
Af$FootpLargeBuildings200m <- rep(0, n)
Af$HasRun <- rep(0, n)
}


#### identify which cities have run, for when it stops
{
    Dr <- dir("HasRun/")
    for(u in 1:length(Dr)){
    v <- as.numeric(gsub(pattern = ".RData",
                    replacement = "",
                    Dr[u])) - 100000
    Af$HasRun[Af$agglosID %in% v] <- 1
    }
}

#### identify city, load buildings as S
#### run metrics and save results
{
    while(min(Af$HasRun==0)){ ### there are missing cities
        k  <- which(Af$HasRun==0)[1]
        timest <- Sys.time()
        #### this files are the buildings dataset for each city
        #### the input is from Google 
    name = paste(Af$agglosID[k], ".csv", sep = "")
    S <- tryCatch(read.csv(name,  header = F),
                  error=function(e){},
                  warning=function(e){},
                  finally = function(e){},
                  silent = T)
    if(is.null(dim(S))){ #### no buildings for that city
        #### output some user info
        {
            cat(" - - - - - - - - - - - - - - - - - - - - \n")
            cat("City - ", k, " - ", as.character(Af$Agglomeration_Name[k]), "\n")
            cat("No buildings for that city \n")
            print(Sys.time() - timest)
            cat(" - - - - - - - - - - - - - - - - - - - - \n")
            Af$HasRun[k] <- 1
            nm <- paste("HasRun/", 100000 + Af$agglosID[k], ".RData", sep = "")
            edf <- data.frame(x = 1)
            save(edf, file = nm)
        }
        S <- data.frame(V1 = c(), V2 = c(), V3 = c(), V4 = c())
        }
    if(dim(S)[1]>0){ ### we found buildings in city
        bw = 0.0025/cos(mean(S$y)/90*pi/2)
        
        #### building metrics
        {
        Af$NBuildings[k] <- dim(S)[1]
        Af$TotalFootprint[k] <- sum(S$area)
        Af$AreaConfidence <- sum(S$area * S$confidence)/sum(S$area)
        Af$MeanBuildingArea[k] <- mean(S$area)
        Af$MedianBuildingArea[k] <- median(S$area)
        Af$VarBuildingArea[k] <- var(S$area)
        Af$AreaBuilt[k] <- sum(S$area)/(Af$Built_up[k]*1000000)
        }
        
        #### polycentrism metrics
        {
        Pol <- PolycentricI(S, bw = bw)
        Af$NBranches[k] <- Pol$NBranches
        Af$Polycentrism[k] <- Pol$Polycentrism
        Af$HerfindahlIndex[k] <- Pol$HerfindahlIndex
        CTree <- Pol$netT
        save(CTree, 
             file = paste("polTrees/PolTree_",
                          Af$agglosID[k],
                          ".RData", sep = ""))
        }
        
        #### fragmentation metrics
        {
        Fr <- FragmentationI(S)
        Af$MeanDistanceBuildings[k] <- Fr$MeanDistanceBuildings
        Af$WMeanDistanceBuildings[k] <- Fr$WMeanDistanceBuildings
        Af$MaxBuildingsDistance[k] <- Fr$MaxBuildingsDistance
        Af$MeanDistanceBuildingsY[k] <- Fr$MeanDistanceBuildingsY
        Af$WMeanDistanceBuildingsY[k] <- Fr$WMeanDistanceBuildingsY
        Af$MaxBuildingsDistanceY[k] <- Fr$MaxBuildingsDistanceY
        Af$MeanDistanceBuildingsG[k] <- Fr$MeanDistanceBuildingsG
        Af$WMeanDistanceBuildingsG[k] <- Fr$WMeanDistanceBuildingsG
        Af$MaxBuildingsDistanceG[k] <- Fr$MaxBuildingsDistanceG
        }
        
        #### density metrics
        {
        Den <- DensityI(S, bw = bw)
        Af$CentreLon[k] <- Den$CentreLon
        Af$CentreLat[k] <- Den$CentreLat
        Af$CentreDensity[k] <- Den$CentreDensity 
        Af$NBuildingsCentre1km[k] <- Den$NBuildingsCentre1km 
        Af$TotalFootprintCentre1km[k] <- Den$TotalFootprintCentre1km
        Af$MeanAreaBuildingsCentre1km[k] <- Den$MeanAreaBuildingsCentre1km
        Af$NBuildingsCentre3km[k] <- Den$NBuildingsCentre3km
        Af$TotalFootprintCentre3km[k] <- Den$TotalFootprintCentre3km
        Af$MeanAreaBuildingsCentre3km[k] <- Den$MeanAreaBuildingsCentre3km
        Af$MeanArea0100m[k] <- Den$MeanArea0100m
        Af$MeanArea0250m[k] <- Den$MeanArea0250m
        Af$MeanArea0500m[k] <- Den$MeanArea0500m
        Af$MeanArea0750m[k] <- Den$MeanArea0750m
        Af$MeanArea1000m[k] <- Den$MeanArea1000m
        Af$MeanArea1500m[k] <- Den$MeanArea1500m
        Af$MeanArea2000m[k] <- Den$MeanArea2000m
        Af$MeanArea2500m[k] <- Den$MeanArea2500m
        Af$MeanArea3000m[k] <- Den$MeanArea3000m
        Af$MeanArea4000m[k] <- Den$MeanArea4000m
        DenIm <- Den$DensityImage
        save(DenIm, 
             file = paste("DensityImages/DensityImage",
                          Af$agglosID[k],
                          ".RData", sep = ""))
        }
        
        ### size distribution
        {
            ADD <- AreaDistDecreasing(S)   
            Af$FootpLargeBuildings050m[k] <- ADD$aBuildings[25]/sum(S$area)
            Af$FootpLargeBuildings100m[k] <- ADD$aBuildings[20]/sum(S$area)
            Af$FootpLargeBuildings200m[k] <- ADD$aBuildings[16]/sum(S$area)
            save(ADD,
                 file = paste("BuildingSize/BSize",
                              Af$agglosID[k],
                              ".RData", sep = ""))
        }
        
        #### save Af files (one per city)
        {
            FAf <- Af[k, ]
            nm <- paste("AfOutput/Africapolis_", Af$agglosID[k], ".csv", sep = "")
            write.csv(FAf, 
                      file = nm,
                      row.names = F)
            nm <- paste("AfOutput/Africapolis_", Af$agglosID[k], ".RData", sep = "")
            save(FAf, file = nmx)
        }
    
        #### create two pdf figures
        {
        nm <- paste("figures/Polycentrism_AgglosID",
                    Af$agglosID[k],
                    ".pdf", sep = "")
        pdf(nm, width = 4, height = 4)
        par(mar = c(.1,.1,.1,.1))
        CityEdge <- owin(poly=list(x=c(min(S$x), max(S$x), max(S$x), min(S$x)),
                               y=c(min(S$y), min(S$y), max(S$y), max(S$y))))
        Den$DensityImage$v <- pmin(Den$DensityImage$v, 4000000000)
        Den$DensityImage$v[1,1] <- 4000000000
        plot(CityEdge, main = "")
        plot(Den$DensityImage, add = T, 
         main = ,
         col = colorRampPalette(c("gray70", "cyan3","gold", "tomato", "tomato3"))(50))
        contour(Den$DensityImage, 
                levels= exp(c(20, 21, 21.5, 22,22.5, 23, 25 )), 
                add = T, 
                col = "white", 
                lwd = .1, 
                drawlabels = F)
        dev.off()
    
    
        nm <- paste("figures/VolumeTree_AgglosID",
                            Af$agglosID[k],
                    ".pdf", sep = "")
        pdf(nm, width = 4, height = 4)
        par(mar = c(.1,.1,.1,.1))
        plotvolu(Pol$netT, 
                 ylab = "", yaxt = "n",
                 ptext=0.002,
                 cex = 1,
                 modelabel=F,symbo="L",
                 paletti = colorRampPalette(c("cyan3"))(50000000),
                 col = "black",
                 colothre = 1,
                 colo=T)
                 
        dev.off()
        }
        
        #### output some user info
        {
        cat(" - - - - - - - - - - - - - - - - - - - - \n")
        cat("City - ", k, " - ", as.character(Af$Agglomeration_Name[k]), "\n")
        cat("number of buildings -- ", dim(S)[1], "\n")
        print(Sys.time() - timest)
        cat(" - - - - - - - - - - - - - - - - - - - - \n")
        Af$HasRun[k] <- 1
        nm <- paste("HasRun/", 100000 + Af$agglosID[k], ".RData", sep = "")
        save(data.frame(x = 1), file = nm)
        }
    }
    }
}

