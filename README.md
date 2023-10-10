# UrbanFragmentation
The repository contains data and code related to African cities' urban metrics and fragmentation.
The data is divided into two separate tables:

## DBMetrics_V2
The table contains the result of the analysis's urban form indicators for all cities. 
The variable to connect to other tables is agglosID

## BuildingSize
The table contains the number of buildings and the total area constructed in a city with buildings smaller than some X. The threshold X is varied with a wide number of values, identified by the names of the columns.

## FragmentationCities GIT
The repository also contains an R Code that reproduces the analysis. It takes as input the Africapolis dataset and the identified buildings within each polygon.
