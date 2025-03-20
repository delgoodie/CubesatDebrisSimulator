// CubesatDebrisModel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Simulator.h"
#include "DataManager.h"

int main()
{
    DataManager dataGenerator;

    Simulator simulator;
    simulator.SetDebris(std::move(dataGenerator.GetData()));
    simulator.Run();
    std::cout << "Debris Detected: " << simulator.DetectionCount << std::endl;
}