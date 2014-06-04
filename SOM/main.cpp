/*
 SOM implementation.
 
 Created by Isidro on 5/25/14.
 Copyright (c) 2014 Isidro Cortes Ciriano. All rights reserved.
 */


#include <iostream>
#include "SOM.h"
#include </usr/local/include/eigen3/Eigen/eigen>
#include <sstream>
#include <fstream>
#include <string>
#include "ClassesSOM.h"

using namespace std;

int main()
{
SOM test;
test.ReadData("/Users/icortes/Desktop/SOM/SOM/SOM/example.csv", 3, 3);
test.InitializeMap(4, 4, 3);
test.Train(10);
test.PrintToCSVFileRowWise("/Users/icortes/Desktop/SOM/SOM/SOM/results.csv", test.SOMMap, 4, 4,3);
return 0;
}

