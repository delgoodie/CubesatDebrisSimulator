#pragma once

#include "CapUtil.h"

class Cubesat
{
public:
	double FieldOfView; // degrees
	double DetectionRange; // km
	
	CoordKep coord;

	Cubesat();
};

