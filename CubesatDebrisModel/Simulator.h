#pragma once
#include "Cubesat.h"

#include <optional>


class Simulator
{
private:
	Cubesat cubesat;
	DebrisList debrisList;
	
	double GravitationalParameter; // (mu) km3/s2
	double TimeStep; // s
	double Duration; // s


public:
	int DetectionCount;

public:
	Simulator();

public:
	void SetDebris(DebrisList InDebrisList) { debrisList = std::move(InDebrisList); }

	void CullDebrisByMinDistance();
	void CullDebrisByIntersection();
	void Run();

private:
	std::optional<double> FindDebrisCubesatPairIntersectionTime(const Debris& debris);
	void SimulateDebrisCubesatPair(const Debris& debris, double StartTime);

	std::vector<double> GenerateTemperatureField(double TimeStep);

public:
	DebrisList& GetDebrisList() {
		return debrisList;
	}
};

