#pragma once
#include "Cubesat.h"




class Simulator
{
private:
	Cubesat cubesat;
	DebrisList debrisList;
	
	float GravitationalParameter; // (mu) km3/s2
	float TimeStep; // s
	float Duration; // s




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

