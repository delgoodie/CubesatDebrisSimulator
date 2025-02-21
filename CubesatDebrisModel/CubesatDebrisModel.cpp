// CubesatDebrisModel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Simulator.h"
#include "DataManager.h"

int main()
{
    DataManager dataGenerator;
    DebrisList debrisList = dataGenerator.GetData();

    Simulator simulator;
    simulator.SetDebris(debrisList);
    simulator.Run();
    std::cout << "Debris Detected: " << simulator.DetectionCount << std::endl;
}