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
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using std::chrono::seconds;



int main()
{
    auto t1 = high_resolution_clock::now();

    DataManager dataManager;

    const char* binDataFileName = "Raw Data/DebrisData2023.bin";
    const char* CSVDataFileName = "Raw Data/OutData.csv";
    const char* LeoDataFileName = "Raw Data/debris_data.json";
    const char* MocatDataFileName = "Raw Data/2023.csv";
    const char* RScriptCsvData = "Raw Data/simulated_debris1.csv";

    const char* culledCsvDataFileName = "Raw Data/CulledDataRdata.csv";
    const char* culledBinDataFileName = "Raw Data/CulledDataRdata.bin";


    // DebrisList debrisList = dataManager.FetchData_bin(culledBinDataFileName);
    // DebrisList debrisList = dataManager.FetchData_bin(culledBinDataFileName);
    // DebrisList debrisList = dataManager.FetchData_bin(culledBinDataFileName);
    // DebrisList debrisList = dataManager.FetchData_csv(culledCsvDataFileName);
    DebrisList debrisList = dataManager.FetchData_csv(RScriptCsvData);

    // DebrisList debrisList = dataManager.GenData_Leo(120e6);
    // dataManager.WriteData_bin(binDataFileName, debrisList);

    // const char* R_script_generated_data = "Raw Data/RScriptData.csv";
    // DebrisList debrisList = dataManager.FetchData_csv(R_script_generated_data);
    // DebrisList debrisList = dataManager.FetchData_csv(culledCsvDataFileName);
    // dataManager.WriteData_bin(culledBinDataFileName , debrisList);


    Simulator simulator;
    simulator.SetDebris(std::move(debrisList));
    simulator.CullDebrisByMinDistance();
    dataManager.WriteData_bin(culledBinDataFileName, simulator.GetDebrisList());
    dataManager.WriteData_csv(culledCsvDataFileName, simulator.GetDebrisList());

    simulator.Run();

    auto t2 = high_resolution_clock::now();
    auto s_int = duration_cast<seconds>(t2 - t1);

    std::cout << "Total application runtime: " << s_int.count() << "s" << std::endl;

    return 0;
}
