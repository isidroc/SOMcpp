//
//  Classes_SOM.h
//  SOM
//
//  Created by Isidro on 6/4/14.
//  Copyright (c) 2014 Isidro Cortes Ciriano. All rights reserved.
//

#ifndef SOM_Classes_SOM_h
#define SOM_Classes_SOM_h

#include </usr/local/include/eigen3/Eigen/eigen>
#include "SOM.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

//typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrix2D;
typedef Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> MatrixOfVectors;


//////////////////////////////////////////////
// Euclidean distance
//////////////////////////////////////////////
class EuclideanDistance {
public:
    float EvaluateKernel(Eigen::VectorXd v1, Eigen::VectorXd v2){
        float KernelValue = (v1 - v2).norm();
        return KernelValue;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////




 //////////////////////////////////////////////
 // Radial Kernel
 //////////////////////////////////////////////
 class RadialKernel {
 public:
 float sigma;
 float EvaluateKernel(Eigen::VectorXd v1, Eigen::VectorXd v2, float sigma){
 float KernelValue = exp(-((v1 - v2).norm() / sigma));
 return KernelValue;
 }
 };
 ////////////////////////////////////////////////////////////////////////////////////////////
 

 //////////////////////////////////////////////
 // Read Input Data
 //////////////////////////////////////////////
class GetInputVectors {
public:
     
 bool fileExist(const char *fileName){
    std::ifstream infile(fileName);
    return infile.good();
 };
     
 float MaxValueInputData;
 MatrixOfVectors  InputVectors;

 void ReadData(const char *filename, int NumberVectors, int InputVectorSize){
 MatrixOfVectors InputVectors(NumberVectors, 1);
 fileExist(filename);
 std::ifstream data(filename);
 std::string line;
 if (data.is_open()){std::cout << "Opening input file..\n\n" << std::endl;}
 else{printf("Cannot open input file.. exiting..\n");}
 
 int rowIndex = 0;
 float MaxValueInputData=0;
 
 while(std::getline(data,line)){
     int colIndex = 0;
     std::stringstream  lineStream(line);
     std::string  cell;
     Eigen::VectorXd InputVector(InputVectorSize);
 
        while(std::getline(lineStream,cell,',')){
            InputVector(colIndex) =std::stof(cell);
            if (std::stof(cell) > MaxValueInputData){MaxValueInputData = std::stof(cell);}
            colIndex++;
        }
     InputVectors(rowIndex,0) = InputVector;
     rowIndex++;
 }
     GetInputVectors::InputVectors = InputVectors;
 }
 };
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
 
 
 //////////////////////////////////////////////
 // SOM Map " Carte "
 //////////////////////////////////////////////
 class Initialize2DMap{
 public:
 MatrixOfVectors SOMMap;
 float MaxValueInputData; // the maximum value in the input data. To have the same range of values in the SOM
 int xsize, ysize;
     
 void InitializeMap(int xsize, int ysize, int InputVectorSize){
     Initialize2DMap::xsize = xsize;
     Initialize2DMap::ysize = ysize;

 //MaxValueInputData = MaxValueInputData;
 MatrixOfVectors SOMMap(xsize,ysize);
     for (int i=0; i<xsize; i++){
        for (int j=0; j<ysize; j++){
            Eigen::VectorXd now = Eigen::VectorXd::Random(InputVectorSize);
            now = now * MaxValueInputData;
            SOMMap(i,j) = now;
        }
     }
     Initialize2DMap::SOMMap = SOMMap;
     Initialize2DMap::xsize = xsize;
     Initialize2DMap::ysize = ysize;
 }
};
 ///////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////
// SOM Learning Function
//////////////////////////////////////////////
class GaussianDistribution {
public:
    float Gaussian(float x, float sigma){
        float value = exp(-(x/sigma));
        return value;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////
// SOM Learning Function
//////////////////////////////////////////////
// Learning (or updating) Function
//      BMU(t+1) = BMU(t) + NeighbourhoodFunction(t)*[InputVector - BMU(t)]
// NeighbourhoodFunction
//      NeighbourhoodFunction(t) = alpha(t) + RadialKernel(InputVector, BMU(t)), where sigma equals alpha(t)^2
// The alpha values at eact t value will be updated with a Gaussian function, with sigma value will be equal to 1/5*(NIters)
// So alpha(t) = exp( - t/5)

class SOMNeighbourhoodFunction : public RadialKernel {
public:
    float NeighbourhoodFunction(float alpha, Eigen::VectorXd input, Eigen::VectorXd bmu){
        float value = alpha + EvaluateKernel(input, bmu, pow(alpha,2));
        return value;
    }
};

///////

class SOMLearningFunction : public SOMNeighbourhoodFunction {
public:
    Eigen::VectorXd LearnFunction(Eigen::VectorXd input, Eigen::VectorXd bmu, float alpha){  //bmu +
        Eigen::VectorXd value =  bmu + NeighbourhoodFunction(alpha, input, bmu)*(input-bmu);
        return value;
    }
};

////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////
/// Print Learing Results
//////////////////////////////////////////////
class PrintSOM {
    public :
    void PrintToCSVFileRowWise (const char *filename, MatrixOfVectors SOMMap, int xsize, int ysize, int InputVectorSize){
        
        std::ofstream file(filename);
        if (file.is_open())
        {
            for (int i=0; i<xsize; i++){  // xsize is considered the rows of the SOM (i or RowMap below), and ysize the columns (j or ColMap below)
                for (int j=0; j<ysize; j++){
                    for (int k=0; k<InputVectorSize; k++){
                        file << SOMMap(i,j)(k) << ",";
                    }
                file << '\n';
                }
            }
        }
    }
};

////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////
/// Train the SOM
//////////////////////////////////////////////
 
class SOMTrain : public GetInputVectors, public Initialize2DMap, public GaussianDistribution,
public SOMLearningFunction, public PrintSOM{
 public:
 int Iters;
 int BestXMap, BestYMap = 0;
 float Score, ScoreNow = 0;
 float alpha;
     
 void Train(int Iters){  //, MatrixOfVectors InputVectors, int xsize, int ysize, MatrixOfVectors SOMMap
     // Initialize alpha
     SOMTrain::alpha = 5.0;
     
     for (int i=0; i<Iters; i++){
         for (int j=0; j<SOMTrain::InputVectors.rows();j++){ //XX the order should be random
             // Find the BMU
             for (int RowMap=0;RowMap<SOMTrain::xsize;RowMap++){
                 for (int ColMap=0;ColMap<ysize;ColMap++){
                     ScoreNow = EvaluateKernel(SOMMap(RowMap,ColMap), InputVectors(j,0),1.0);
                     if (SOMTrain::ScoreNow > Score){SOMTrain::Score = ScoreNow; SOMTrain::BestXMap=RowMap, SOMTrain::BestYMap=ColMap;}
                  }
             }
             // Update the BMUs
             for (int RowMap=0;RowMap<SOMTrain::xsize;RowMap++){
                 for (int ColMap=0;ColMap<ysize;ColMap++){
                     Eigen::VectorXd BMUUpdated = LearnFunction(InputVectors(j,0), SOMMap(RowMap,ColMap), SOMTrain::alpha);
                     
                 }
             }
             
         }
     // Update alpha
         SOMTrain::alpha = Gaussian(i,5);
         //std::cout << alpha << SOMTrain::alpha << std::endl;
 }
}
};
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
