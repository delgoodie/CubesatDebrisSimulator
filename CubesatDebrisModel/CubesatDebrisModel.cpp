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


float frand(float LO, float HI) {
    return  LO + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX / (HI - LO)));
}

int main()
{
    // DataManager dataGenerator;

    float OverallMin = 0;
    CoordKep BestCoordA, BestCoordB;

    std::vector<glm::vec3> SamplesA;
    std::vector<glm::vec3> SamplesB;

    for (int i = 0; i < 10000; i++) 
    {
        SamplesA.clear();
        SamplesB.clear();

        CoordKep CoordA = { 8000.f, 0.f, CapUtil::Deg2Rad(frand(0.f, 90.f)), CapUtil::Deg2Rad(frand(0.f, 90.f)), CapUtil::Deg2Rad(frand(0.f, 90.f)), 0.f };
        CoordKep CoordB = { 8000.f, 0.f, CapUtil::Deg2Rad(frand(0.f, 90.f)), CapUtil::Deg2Rad(frand(0.f, 90.f)), CapUtil::Deg2Rad(frand(0.f, 90.f)), 0.f };

        float MinDist = CapUtil::MinDistance(
            CoordA, CoordB,
            SamplesA, 
            SamplesB
        );
        if (MinDist > OverallMin) {
            OverallMin = MinDist;
            BestCoordA = CoordA;
            BestCoordB = CoordB;
        }
    }



    std::cout << OverallMin << std::endl;
    std::cout << CapUtil::Rad2Deg(BestCoordA.i) << "  " << CapUtil::Rad2Deg(BestCoordA.o) << "  " << CapUtil::Rad2Deg(BestCoordA.w) << std::endl;
    std::cout << CapUtil::Rad2Deg(BestCoordB.i) << "  " << CapUtil::Rad2Deg(BestCoordB.o) << "  " << CapUtil::Rad2Deg(BestCoordB.w) << std::endl;
    SamplesA.clear();
    SamplesB.clear();
    float MinDist2 = CapUtil::MinDistance(
        BestCoordA, BestCoordB,
        SamplesA,
        SamplesB
    );


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


    // Simulator simulator;
    // simulator.SetDebris(std::move(dataGenerator.GetData()));
    // simulator.Run();
    // std::cout << "Debris Detected: " << simulator.DetectionCount << std::endl;

    return 0;
}