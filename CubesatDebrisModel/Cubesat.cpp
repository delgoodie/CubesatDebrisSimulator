#include "Cubesat.h"

Cubesat::Cubesat()
{
	FieldOfView = 60.f;
	DetectionRange = 1;

	coord = { 6928, 0, CapUtil::Deg2Rad(37.556), CapUtil::Deg2Rad(200), 0, 0 };
	// CapUtil::CC_to_CK({ glm::vec3(2000.f, 0.f, 0.f), glm::vec3(0.f, 7.5f, 0.f) });
}
