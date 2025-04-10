#include "Cubesat.h"

Cubesat::Cubesat()
{
	FieldOfView = 60.;
	DetectionRange = 1;

	coord = { 6928, 0, CapUtil::Deg2Rad(37.556), CapUtil::Deg2Rad(200), 0, 0 };
	// CapUtil::CC_to_CK({ vec3(2000., 0., 0.), vec3(0., 7.5, 0.) });
}
