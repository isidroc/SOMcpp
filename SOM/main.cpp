/*
 SOM implementation.
 
 Created by Isidro on 5/25/14.
 Copyright (c) 2014 Isidro Cortes Ciriano. All rights reserved.
 */



#include </usr/local/include/eigen3/Eigen/eigen>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include "SOM.h"
#include "ClassesSOM.h"

using namespace std;

int main()
{
SOM test;
test.ReadData("/Users/icortes/Desktop/SOMcpp/SOM/colorsInput.csv", 3000, 3);
test.InitializeMap(250, 250);
test.Train(25000 + 250000);
test.PrintToCSVFileRowWise("/Users/icortes/Desktop/SOMcpp/SOM/ColorsResults.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);
    
    //std::cout << test.MaxValueInputData;
    return 0;
}

