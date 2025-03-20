#include "Cubesat.h"

Cubesat::Cubesat()
{
	FieldOfView = 60.f;
	DetectionRange = 1;

	coord = CapUtil::CC_to_CK({ glm::vec3(2000.f, 0.f, 0.f), glm::vec3(0.f, 7.5f, 0.f) });
}
