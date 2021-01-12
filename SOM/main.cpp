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
#include <chrono>

using namespace std;

int main()
{
SOM test;
test.ReadData("C:/Users/4ever/Documents/SOMcpp-master/SOM/example.csv", 13, 3);
//test.ReadData("C:/Users/4ever/Documents/GitHub/SOMcpp/SOM/ThreeColors.txt", 3, 3);
test.InitializeMap(100, 100);
test.PrintToCSVFileRowWise("C:/Users/4ever/Documents/GitHub/SOMcpp/SOM/SOMColorsBeforeTraining2.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);

    //test.SigmaNeighbouringInitial = 100;
    //test.SigmaNeighbourhoodFinal = 1;
    //test.LearningRateInitial = 1.0;
    //test.LearningRateFinal = 1.0;
std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
test.Train(1000);
std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
std::cout << "Time ellapsed = " << (std::chrono::duration_cast<std::chrono::nanoseconds> (end - begin).count())/1000000000.0 << "[s]" << std::endl;
test.PrintToCSVFileRowWise("C:/Users/4ever/Documents/GitHub/SOMcpp/SOM/ColorsResults2.csv", test.SOMMap, test.xsize, test.ysize,test.InputVectorSize);
    
    //std::cout << test.MaxValueInputData;
    
    //std::cout <<  "L :" << test.LearningRate << "\n";
    //std::cout <<  "Sigma Neighbours:" << test.SigmaNeighbouring << "\n";
    return 0;
}

