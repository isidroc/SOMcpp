/*
 SOM implementation.
 
 Created by Isidro on 5/25/14.
 Copyright (c) 2014 Isidro Cortes Ciriano. All rights reserved.
 */





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
//test.ReadData("C:/Users/4ever/Documents/SOMcpp-master/SOM/example.csv", 13, 3);
test.ReadData("C:/Users/4ever/Documents/SOMcpp-master/SOM/ThreeColors.txt", 3, 3);
test.InitializeMap(50, 50);
test.PrintToCSVFileRowWise("C:/Users/4ever/Documents/SOMcpp-master/SOM/SOMColorsBeforeTraining2.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);

    test.SigmaNeighbouringInitial = 50;
    test.SigmaNeighbourhoodFinal = 1;
    test.LearningRateInitial = 0.6;
    test.LearningRateFinal = 0.25;
    
test.Train(200);
test.PrintToCSVFileRowWise("C:/Users/4ever/Documents/SOMcpp-master/SOM/ColorsResults2.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);
    
    //std::cout << test.MaxValueInputData;
    
    std::cout << "lambda :" << test.lambda << "\n";
    std::cout <<  "L0 :" << test.LearningRateInitial << "\n";
    std::cout <<  "Sigma0 Neighbours:" << test.SigmaNeighbouringInitial << "\n";
    //std::cout <<  "L :" << test.LearningRate << "\n";
    //std::cout <<  "Sigma Neighbours:" << test.SigmaNeighbouring << "\n";
    return 0;
}

