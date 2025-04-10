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
    DataManager dataManager;

    const char* binDataFileName = "Raw Data/DebrisData2023.bin";
    const char* CSVDataFileName = "Raw Data/OutData.csv";
    const char* LeoDataFileName = "Raw Data/debris_data.json";
    const char* MocatDataFileName = "Raw Data/2023.csv";

    const char* culledCsvDataFileName = "Raw Data/CulledData.csv";
    const char* culledBinDataFileName = "Raw Data/CulledData.bin";


    // DebrisList debrisList = FetchData_bin(binDataFileName);

    //DebrisList debrisList = dataManager.GenData_Leo(120e6);
    // dataManager.WriteData_bin(binDataFileName, debrisList);


    DebrisList debrisList = dataManager.FetchData_bin(culledBinDataFileName);
    // DebrisList debrisList = dataManager.FetchData_csv(culledCsvDataFileName);
    // dataManager.WriteData_bin(culledBinDataFileName , debrisList);


    Simulator simulator;
    simulator.SetDebris(std::move(debrisList));
    // simulator.CullDebrisByMinDistance();
    // dataManager.WriteData_bin(culledBinDataFileName, simulator.GetDebrisList());
    // dataManager.WriteData_csv(culledCsvDataFileName, simulator.GetDebrisList());


    // if (!dataManager.HasDataFile(binDataFileName))
    // {
    // }

    simulator.Run();
    std::cout << "Debris Detected: " << simulator.DetectionCount << std::endl;

    return 0;
}
