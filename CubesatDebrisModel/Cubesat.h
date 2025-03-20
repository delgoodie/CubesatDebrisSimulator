#pragma once

#include "CapUtil.h"

class Cubesat
{
public:
	float FieldOfView; // degrees
	float DetectionRange; // km
	
	CoordKep coord;

	Cubesat();
};

