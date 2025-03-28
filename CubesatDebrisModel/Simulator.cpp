#include "Simulator.h"

#include <iostream>
#include "CapUtil.h"
#include <vector>

#define PI 3.1415926535f

void UpdateOrbitalBody(CoordCar& coord, float GravitationalParameter, float TimeStep)
{
    glm::vec3 PrevPosition = coord.pos;
    
    float Radius = glm::length(coord.pos); // sqrt(Position[0] * Position[0] + Position[1] * Position[1] + Position[2] * Position[2]);

    glm::vec3 Acceleration = coord.pos * -GravitationalParameter / (Radius * Radius * Radius);
    coord.vel += Acceleration * TimeStep;
    coord.pos += coord.vel *TimeStep;

    // std::cout << PrevPosition[0] << ", " << PrevPosition[1] << ", " << PrevPosition[2] << " --> " << Position[0] << ", " << Position[1] << ", " << Position[2] << std::endl;
}


void UpdateOrbitalBody(CoordKep& coord, float GravitationalParameter, float TimeStep)
{
    coord.m += fmod(sqrt(GravitationalParameter / (coord.a * coord.a * coord.a)) * TimeStep, 2 * PI);
}



Simulator::Simulator()
{
	GravitationalParameter = 398600.4418f;
    TimeStep = 1.f; // 1.f / 60.f;
	Duration = 10 * 1;

    DetectionCount = 0;
}




void Simulator::CullDebrisByMinDistance()
{
    std::cout << "Culling debris by distance" << std::endl;

    int InitialSize = (int)debrisList.num;

    std::vector<Debris> DebrisInRange;
    for (int i = 0; i < debrisList.num; i++) 
    {
        if ((i % int(debrisList.num / 100)) == 0)
        {
            std::cout << i << "..." << std::endl;
        }
        float MinDist = CapUtil::MinDistanceBetweenEllipses(cubesat.coord, debrisList[i].coord);
        if (MinDist <= cubesat.DetectionRange * 1.5f)
        {
            DebrisInRange.push_back(debrisList[i]);
        }
    }

    debrisList = DebrisList(DebrisInRange.size());
    for (int i = 0; i < DebrisInRange.size(); i++)
    {
        debrisList[i] = DebrisInRange[i];
    }

    int FinalSize = (int)debrisList.num;

    std::cout << "Culled " << InitialSize - FinalSize << " debris out of " << InitialSize << " initial debris (" << (1.f - ((float)FinalSize / float(InitialSize))) * 100.f << "%)" << std::endl;
}

void Simulator::Run()
{
    CullDebrisByMinDistance();


    if (debrisList.num == 0) 
    {
        std::cout << "Can't Run Simulation: Debris Position / Velocity unset" << std::endl;
        return;
    }

    DetectionCount = 0;
    int Steps = int(Duration / TimeStep);

    std::cout << "Start Simulation (steps=" << Steps << ")" << std::endl;

    for (int i = 0; i < Steps; i++) 
    {

        UpdateOrbitalBody(cubesat.coord, GravitationalParameter, TimeStep);
        glm::vec3 cubesatPos = CapUtil::CK_to_CC(cubesat.coord).pos;

        
        for (int d = 0; d < debrisList.num; d++)
        {
            // if mod(i, floor(this.num_debris / 10)) == 0
            //    fprintf("Debris %d\n", i)

            UpdateOrbitalBody(debrisList[d].coord, GravitationalParameter, TimeStep);

            glm::vec3 debrisPos = CapUtil::CK_to_CC(debrisList[i].coord).pos;

            glm::vec3 difference = cubesatPos - debrisPos;

            float DistanceSquared = glm::dot(difference, difference);
            // std::cout << DistanceSquared << std::endl;
            bool bDetected = DistanceSquared < (cubesat.DetectionRange * cubesat.DetectionRange);

            DetectionCount += bDetected;
        }
        std::cout << "Simulator Iteration" << std::endl;
    }
    std::cout << "Simulator Finished" << std::endl;
}
