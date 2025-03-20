#include "DataManager.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <filesystem>
#include <json/json.hpp>


using json = nlohmann::json;


DebrisList DataManager::GetData()
{
    const bool bUseCachedFile = false;
    const bool bWriteToCSV = true;

    const char* binDataFileName = "Raw Data/DebrisData2023.bin";
    const char* CSVDataFileName = "Raw Data/OutData.csv";
    if (!bUseCachedFile || !std::filesystem::exists(binDataFileName))
    {
        // DebrisList debrisList = FetchMocatData_csv("Raw Data/2023.csv");
        DebrisList debrisList = FetchLeoData_json("Raw Data/debris_data.json");

        if (bUseCachedFile)
        {
            WriteDataToBinFile(binDataFileName, debrisList);
        }

        if (bWriteToCSV) 
        {
            WriteDataToCSV(CSVDataFileName, debrisList);
        }


        return debrisList;
    }

    DebrisList debrisList = FetchData_bin(binDataFileName);

    return debrisList;
}


void DataManager::WriteDataToBinFile(const char* FileName, DebrisList& debrisList)
{
    std::ofstream outFile(FileName, std::ios::binary);
    if (!outFile) {
        std::cerr << "Error opening file for writing!\n";
        return;
    }
    outFile.write(reinterpret_cast<const char*>(debrisList.GetRawData()), debrisList.Num() * sizeof(Debris));
    outFile.close();
}

void DataManager::WriteDataToCSV(const char* FileName, DebrisList& debrisList)
{
    std::ofstream outFile(FileName);

    if (!outFile) {
        std::cerr << "Error opening file for writing!\n";
        return;
    }
    
    outFile << "a, e, i, o, w, m, rcs\n";

    for (int i = 0; i < debrisList.Num(); i++) 
    {
        outFile << debrisList[i].coord.a << ", " << debrisList[i].coord.e << ", " << debrisList[i].coord.i << ", " << debrisList[i].coord.o << ", " << debrisList[i].coord.w << ", " << debrisList[i].coord.m << ", " << debrisList[i].rcs << "\n";
    }

    outFile.close();
    std::cout << "CSV file written successfully!\n";
}

DebrisList DataManager::FetchMocatData_csv(const char* FileName)
{
    std::ifstream file(FileName);
    if (!file.is_open()) {
        std::cerr << "Error opening file" << std::endl;
        return DebrisList();
    }

    std::vector<Debris> debrisVec;

    std::string line;
    int LineNumber = -1;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string cell;

        LineNumber++;

        if (LineNumber == 0) continue;

        Debris debris;

        int col = 0;
        while (std::getline(ss, cell, ',')) {
            if (col == 12)
            {
                debris.coord.e = std::stof(cell.substr(1, cell.length() - 2));
            }
            else if (col == 13) 
            {
                debris.coord.i = CapUtil::Deg2Rad(std::stof(cell.substr(1, cell.length() - 2)));
            }
            else if (col == 14)
            {
                debris.coord.o = CapUtil::Deg2Rad(std::stof(cell.substr(1, cell.length() - 2)));
            }
            else if (col == 15)
            {
                debris.coord.w = CapUtil::Deg2Rad(std::stof(cell.substr(1, cell.length() - 2)));
            }
            else if (col == 16)
            {
                debris.coord.m = CapUtil::Deg2Rad(std::stof(cell.substr(1, cell.length() - 2)));
            }
            else if (col == 25)
            {
                debris.coord.a = std::stof(cell.substr(1, cell.length() - 2));
            }
            else if (col == 30)
            {
                std::string SizeStr = cell.substr(1, cell.length() - 2);
                bool bIsLarge = SizeStr._Starts_with("L");
                bool bIsMedium = SizeStr._Starts_with("M");
                debris.rcs = bIsLarge ? 1.f : bIsMedium ? .5f : .1f;
                break;
            }
            debris.rcs = 0.f;
            col++;
        }

        debrisVec.push_back(debris);
    }

    file.close();

    DebrisList debrisList(debrisVec.size());

    for (int i = 0; i < debrisVec.size(); i++) 
    {
        debrisList[i] = debrisVec[i];
    }

    return debrisList;
}

DebrisList DataManager::FetchLeoData_json(const char* FileName)
{
    std::cout << "Loading debris from " << FileName << std::endl;

    std::ifstream data_file(FileName, std::ifstream::binary);
    json jsonData = json::parse(data_file);

    int NumDebris = 0; // jsonData["objects"].size();
    int NumTotal = 0;

    std::vector<std::tuple<std::string, int>> Types;
    for (auto orbitObj : jsonData["objects"])
    {
        NumTotal++;
        if (orbitObj["type"] == "debris")
        {
            NumDebris++;
        }
        bool bFoundType = false;
        for (int i = 0; i < Types.size(); i++) {
            if (std::get<std::string>(Types[i]) == orbitObj["type"]) {
                bFoundType = true;
                std::get<int>(Types[i])++;
            }
        }
        if (!bFoundType) {
            Types.push_back({ orbitObj["type"] , 1 });
        }
    }
    for (int k = 0; k < Types.size(); k++)
    {
        std::cout << "Loaded " << std::get<int>(Types[k]) << " " << std::get<std::string>(Types[k]) << std::endl;
    }
    std::cout << NumTotal << "Total Objects" << std::endl;

    DebrisList debrisList(NumDebris);

    int i = 0;
    for (auto orbitObj : jsonData["objects"])
    {
        if (orbitObj["type"] == "debris")
        {
            glm::vec3 pos = { (float)orbitObj["position"][0], (float)orbitObj["position"][1], (float)orbitObj["position"][2] };
            glm::vec3 vel = { (float)orbitObj["velocity"][0], (float)orbitObj["velocity"][1], (float)orbitObj["velocity"][2] };
            debrisList[i].coord = CapUtil::CC_to_CK({ pos, vel });
            debrisList[i].rcs = 0;
            i++;
        }
    }

    return debrisList;
}

DebrisList DataManager::FetchData_bin(const char* FileName)
{
    std::ifstream inFile(FileName, std::ios::binary);
    if (!inFile) {
        std::cerr << "Error opening file for reading!\n";
        return DebrisList(0);
    }

    // Get file size
    inFile.seekg(0, std::ios::end);
    size_t fileSize = inFile.tellg();
    inFile.seekg(0, std::ios::beg);

    DebrisList debrisList(fileSize / sizeof(Debris));

    inFile.read(reinterpret_cast<char*>(debrisList.GetRawData()), fileSize);
    inFile.close();

    return debrisList;
}