//
//  Classes_SOM.h
//  SOM
//
//  Created by Isidro on 6/4/14.
//  Copyright (c) 2014 Isidro Cortes Ciriano. All rights reserved.
//

#ifndef SOM_Classes_SOM_h
#define SOM_Classes_SOM_h

#include <time.h>
#include </eigen/Eigen/eigen>
#include "SOM.h"
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>

//typedef Eigen::Matrix<float, Eigen::Dynamic, Eigen::Dynamic> Matrix2D;
typedef Eigen::Matrix<Eigen::VectorXd, Eigen::Dynamic, Eigen::Dynamic> MatrixOfVectors;

//////////////////////////////////////////////
// Random Generator
//////////////////////////////////////////////
class RandomGenerator {
    public :
    RandomGenerator(){
        srand(time(NULL));
    }
   float GenerateRandomNumber(int RandMax){
       float value = rand() % RandMax;
        return value;
    }
};
//////////////////////////////////////////////


//////////////////////////////////////////////
// Euclidean distance
//////////////////////////////////////////////
class EuclideanDistance {
public:
    float EvaluateDistance(Eigen::VectorXd v1, Eigen::VectorXd v2){
        //float KernelValue = (v1 - v2).norm();
        double distance = 0;
        for (int i=0; i<v2.cols(); ++i)
        {
            distance += (v2(i) - v1(i)) * (v2(i) - v1(i));
        }
        return sqrt(distance);
       // return KernelValue;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////
// Radial Kernel Distance
//////////////////////////////////////////////
class RadialKernelDistance {
public:
    float sigma;
    float EvaluateDistance(Eigen::VectorXd v1, Eigen::VectorXd v2, float sigma){
        float KernelValue = exp(-((v1 - v2).norm() / sigma));
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
    float MaxValueInputData=0;
    MatrixOfVectors  InputVectors;
    int InputVectorSize;
    int NumberInputVectors;
    
 bool fileExist(const char *fileName){
    std::ifstream infile(fileName);
    return infile.good();
 };

 void ReadData(const char *filename, int NumberVectors, int InputVectorSize){
 MatrixOfVectors InputVectors(NumberVectors, 1);
 
 GetInputVectors::InputVectorSize = InputVectorSize;
 GetInputVectors::NumberInputVectors = NumberVectors;
 //float MaxValueInputData;
 fileExist(filename);
 std::ifstream data(filename);
 std::string line;
 if (data.is_open()){std::cout << "Opening input file..\n\n" << std::endl;}
 else{printf("Cannot open input file.. exiting..\n");}
 
 int rowIndex = 0;

 
 while(std::getline(data,line)){
     int colIndex = 0;
     std::stringstream  lineStream(line);
     std::string  cell;
     Eigen::VectorXd InputVector(InputVectorSize);
 
        while(std::getline(lineStream,cell,',')){
            
            InputVector(colIndex) =std::stof(cell);
            //std::cout << cell <<std::endl;
            if (std::stof(cell) > GetInputVectors::MaxValueInputData){GetInputVectors::MaxValueInputData = std::stof(cell);}
            colIndex++;
        }
     InputVectors(rowIndex,0) = InputVector;
     rowIndex++;
 }
     data.close();
     printf("Closing input file..\n");
     GetInputVectors::InputVectors = InputVectors;
     
 }
 };
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
 
 
 //////////////////////////////////////////////
 // SOM Map " Carte "
 //////////////////////////////////////////////
class Initialize2DMap : public GetInputVectors{
 public:
 MatrixOfVectors SOMMap;
 float MaxValueInputData; // the maximum value in the input data. To have the same range of values in the SOM
 int xsize, ysize;
     
 void InitializeMap(int xsize, int ysize){
     Initialize2DMap::xsize = xsize;
     Initialize2DMap::ysize = ysize;
     MatrixOfVectors SOMMap(xsize,ysize);
     
     for (int i=0; i<xsize; i++){
        for (int j=0; j<ysize; j++){
            Eigen::VectorXd now = Eigen::VectorXd::Random(InputVectorSize);
            now  = (now.cwiseAbs() / now.cwiseAbs().maxCoeff()) * GetInputVectors::MaxValueInputData; // random values between 0 and 1
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
        float value = exp(-((float)x/(float)sigma));
        return value;
    }
};
////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////
// SOM Learning Function
//////////////////////////////////////////////
// Learning (or updating) Function
//      BMU(t+1) = BMU(t) + LearningRate(t) * NeighbourhoodFunction(t)*[InputVector - BMU(t)]
// NeighbourhoodFunction
//
class SOMNeighbourhoodFunction  {
public:
    float NeighbourhoodFunction(int BMUX, int BMUY, int BMUNowX, int BMUNowY, float SigmaNeighbouring){
        float value = (float)(pow(BMUX - BMUNowX,2) + pow(BMUY - BMUNowY,2)) / (float)(2*pow(SigmaNeighbouring,2)) ;
        value = exp(-value);
       // value = exp(value);
       // std::cout  << (float)(pow(BMUX - BMUNowX,2) + pow(BMUY - BMUNowY,2)) << (float)2*pow(SigmaNeighbouring,2) << value << "\n";
        //std::cout << value << "\n";
        return value;
    }
};

///////

class SOMLearningFunction : public SOMNeighbourhoodFunction {
public:
    Eigen::VectorXd LearnFunction(Eigen::VectorXd BMUNow, Eigen::VectorXd Input, int BMUX,int BMUY,int BMUNowX, int BMUNowY,float SigmaNeighbouring, float LearningRate){
        Eigen::VectorXd value = BMUNow + LearningRate * NeighbourhoodFunction(BMUX,BMUY,BMUNowX,BMUNowY,SigmaNeighbouring)*(Input-BMUNow);
        return value;
    }
};

///////

class UpdateSOMLearningRate {
public:
    float UpdateLearningRate(float Iter,int maxIters){
        float value = 0.1*exp(-((float)Iter/maxIters));
        return value;
    }
};

///////

class UpdateSOMSigmaNeighbouring {
public:
    float UpdateSigmaNeighbouring(float Iter,int maxIters,int radius){
        float value = radius*exp(-((float)Iter/(maxIters/(log(radius)))));
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
 
class SOMTrain :  public Initialize2DMap, public GaussianDistribution,
public SOMLearningFunction, public PrintSOM, public EuclideanDistance, public UpdateSOMLearningRate, public UpdateSOMSigmaNeighbouring, public RandomGenerator {
 public:
 int Iters;
 int BestXMap, BestYMap = 0;
 float Score, ScoreNow = 0;
    

     
 void Train(int Iters){  //, MatrixOfVectors InputVectors, int xsize, int ysize, MatrixOfVectors SOMMap
     // Initialize alpha
     
     float LearningRate;
     float SigmaNeighbouring;
     //std::cout <<  SOMTrain::SigmaNeighbouring[1] <<std::endl;
     int finned =0;
     Eigen::VectorXd BMU;
     int bmuX;
     int bmuY;
     while(finned<Iters){
         LearningRate = UpdateLearningRate( finned,Iters);
        SigmaNeighbouring = UpdateSigmaNeighbouring(finned,Iters,SOMTrain::xsize);
        int j = GenerateRandomNumber(SOMTrain::NumberInputVectors);
        SOMTrain::Score = 0;
        for (int RowMap=0;RowMap<SOMTrain::xsize;RowMap++){
            for (int ColMap=0;ColMap<SOMTrain::ysize;ColMap++){
                ScoreNow = EvaluateDistance(Initialize2DMap::SOMMap(RowMap,ColMap), InputVectors(j,0)); //XX ,1.0 sigma
                if (ScoreNow > SOMTrain::Score){SOMTrain::Score = ScoreNow; SOMTrain::BestXMap=RowMap; SOMTrain::BestYMap=ColMap;}
            }
        }
        //BMU = Initialize2DMap::SOMMap(SOMTrain::BestXMap,SOMTrain::BestYMap);
        bmuX = SOMTrain::BestXMap;
        bmuY = SOMTrain::BestYMap;
             
        // Update the BMUs
        for (int RowMap=0;RowMap<SOMTrain::xsize;RowMap++){
            for (int ColMap=0;ColMap<ysize;ColMap++){
                
                
                int CoordX = RowMap;
                int CoordY = ColMap;
                // if (CoordX < (SOMTrain::xsize/2)) {CoordX = RowMap;}
            //  if (CoordX >= (SOMTrain::xsize/2)){CoordX = (SOMTrain::xsize/2) - RowMap%(SOMTrain::xsize/2);}
                //if (CoordY < (SOMTrain::ysize/2)) {CoordY = ColMap;}
                //if (CoordY >= (SOMTrain::ysize/2)){CoordY = (SOMTrain::ysize/2) - ColMap%(SOMTrain::xsize/2);}
                //std::cout << CoordX << "\n";// ;<< "  " << CoordY << std::endl;
                
                Eigen::VectorXd BMUUpdated = LearnFunction(Initialize2DMap::SOMMap(RowMap,ColMap),InputVectors(j,0),
                bmuX, bmuY, CoordX, CoordY, SigmaNeighbouring, LearningRate);
                Initialize2DMap::SOMMap(RowMap,ColMap) = BMUUpdated;
                //std::cout << SOMTrain::BestXMap << std::endl;
                
            }
        }
                //j+(i*InputVectors.rows());
                
            //}
        // Update Learning Rate
            
               finned++;
    }
        /*for (i[j]=0; i<Iters; ){
            //for (int j=0; j<SOMTrain::InputVectors.rows();j++){ //XX the order should be random
            
                // Find the BMU
                
        }*/
    }
};
 ///////////////////////////////////////////////////////////////////////////////////////////////////////
#endif
