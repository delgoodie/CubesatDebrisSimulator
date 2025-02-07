// CubesatDebrisModel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Simulator.h"
#include "DebrisDataGenerator.h"

int main()
{
    DebrisDataGenerator dataGenerator;
    dataGenerator.LoadData();

    float* DebrisPosition;
    float* DebrisVelocity;
    int NumDebris;
    dataGenerator.TakeData(DebrisPosition, DebrisVelocity, NumDebris);


    Simulator simulator;
    simulator.SetDebris(DebrisPosition, DebrisVelocity, NumDebris);
    simulator.Run();
    std::cout << "Debris Detected: " << simulator.DetectionCount << std::endl;
}