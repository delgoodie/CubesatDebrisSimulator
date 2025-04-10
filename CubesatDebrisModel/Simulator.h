#pragma once
#include "Cubesat.h"




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
	void Run();


	DebrisList& GetDebrisList() {
		return debrisList;
	}
};

