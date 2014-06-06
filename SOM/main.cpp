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
#include <stdlib.h>     /* srand, rand */


using namespace std;

int main()
{

    SOM test;
test.ReadData("/Users/icortes/Desktop/SOMcpp/SOM/colorsInput.csv", 3000, 3);
//test.ReadData("/Users/icortes/Desktop/SOMcpp/SOM/ThreeColors.txt", 3, 3);
test.InitializeMap(50, 50);
test.PrintToCSVFileRowWise("/Users/icortes/Desktop/SOMcpp/SOM/SOMColorsBeforeTraining.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);

    test.SigmaNeighbouringInitial = 50;
    test.LearningRateInitial = 0.1;
    
test.Train(1000);
test.PrintToCSVFileRowWise("/Users/icortes/Desktop/SOMcpp/SOM/ColorsResults.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);
    
    //std::cout << test.MaxValueInputData;
    
   std::cout <<  test.LearningRate << "\n";
   std::cout <<  test.SigmaNeighbouring << "\n";
    return 0;
}

