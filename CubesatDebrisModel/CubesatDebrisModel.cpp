// CubesatDebrisModel.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "Simulator.h"
#include "DataManager.h"
#include "CapUtil.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>


int main()
{
    // DataManager dataGenerator;

    std::vector<glm::vec3> SamplesA;
    std::vector<glm::vec3> SamplesB;
    float MinDist = CapUtil::MinDistance({ 8000.f, 0.f, CapUtil::Deg2Rad(20.f), CapUtil::Deg2Rad(40.f), 0.f, 0.f }, { 8000.f, 0.f, CapUtil::Deg2Rad(30.f), CapUtil::Deg2Rad(0.f), 0.f, 0.f }, SamplesA, SamplesB);

    std::ofstream outFile("out_test.csv");

    if (!outFile) {
        std::cerr << "Error opening file for writing!\n";
    }

    outFile << "Ax, Ay, Az, Bx, By, Bz\n";

    for (int i = 0; i < SamplesA.size(); i++)
    {
        outFile << SamplesA[i].x << ", " << SamplesA[i].y << ", " << SamplesA[i].z << ", " << SamplesB[i].x << ", " << SamplesB[i].y << ", " << SamplesB[i].z << "\n";
    }

    outFile.close();
    std::cout << "CSV file written successfully!\n";


    std::cout << MinDist << std::endl;

    // Simulator simulator;
    // simulator.SetDebris(std::move(dataGenerator.GetData()));
    // simulator.Run();
    // std::cout << "Debris Detected: " << simulator.DetectionCount << std::endl;

    return 0;
}