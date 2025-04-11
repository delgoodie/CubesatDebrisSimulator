#pragma once

#include "CapUtil.h"

class Cubesat
{
public:
	double FieldOfView; // degrees
	double DetectionRange; // km
	double TransmitPower; // W
	double TransmitAntennaGainFactor; // dbi
	double Frequency; // Hhz
	double Wavelength; // m

	CoordKep coord;

	Cubesat();
};

