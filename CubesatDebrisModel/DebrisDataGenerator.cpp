#include "DebrisDataGenerator.h"

#include "json.hpp"
#include <fstream>
#include <iostream>


using json = nlohmann::json;

DebrisDataGenerator::DebrisDataGenerator()
{
	DebrisPosition = nullptr;
	DebrisVelocity = nullptr;
	NumDebris = 0;
}

void DebrisDataGenerator::FetchData()
{
	
}

void DebrisDataGenerator::CacheData()
{
}

void DebrisDataGenerator::LoadData()
{
	std::cout << "Loading debris from " << "debris_data.json" << std::endl;

	std::ifstream data_file("debris_data.json", std::ifstream::binary);
	json jsonData = json::parse(data_file);

	NumDebris = 0; // jsonData["objects"].size();
	for (auto orbitObj : jsonData["objects"])
	{
		if (orbitObj["type"] == "debris") 
		{
			NumDebris++;
		}
	}
	std::cout << "Loaded " << NumDebris << " Debris" << std::endl;

	DebrisPosition = (float*)malloc(sizeof(float) * NumDebris * 3);
	DebrisVelocity = (float*)malloc(sizeof(float) * NumDebris * 3);

	int i = 0;
	for (auto orbitObj : jsonData["objects"])
	{
		if (orbitObj["type"] == "debris")
		{
			DebrisPosition[i*3+0] = (float)orbitObj["position"][0];
			DebrisPosition[i*3+1] = (float)orbitObj["position"][1];
			DebrisPosition[i*3+2] = (float)orbitObj["position"][2];

			DebrisVelocity[i*3+0] = (float)orbitObj["velocity"][0];
			DebrisVelocity[i*3+1] = (float)orbitObj["velocity"][1];
			DebrisVelocity[i*3+2] = (float)orbitObj["velocity"][2];
			i++;
		}
	}
}

void DebrisDataGenerator::TakeData(float*& OutDebrisPosition, float*& OutDebrisVelocity, int& OutNumDebris)
{
	OutDebrisPosition = DebrisPosition;
	OutDebrisVelocity = DebrisVelocity;
	OutNumDebris = NumDebris;
	DebrisPosition = nullptr;
	DebrisVelocity = nullptr;
	NumDebris = -1;
}
