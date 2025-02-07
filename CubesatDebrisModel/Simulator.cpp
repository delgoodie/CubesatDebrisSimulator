#include "Simulator.h"

#include <iostream>
#include "VectorUtil.h"


#define PI 3.1415926535

void PrintVec(float Vec[3])
{
    std::cout << Vec[0] << ", " << Vec[1] << ", " << Vec[2] << std::endl;
}


void UpdateOrbitalBody(float Position[3], float Velocity[3], float GravitationalParameter, float TimeStep)
{
    float PrevPosition[3];
    VectorUtil::Copy(PrevPosition, Position);
    
    float Radius = sqrt(Position[0] * Position[0] + Position[1] * Position[1] + Position[2] * Position[2]);


    float Acceleration[3];
    VectorUtil::ScalarMultiply(Acceleration, Position, -GravitationalParameter / (Radius * Radius * Radius) * TimeStep);
    VectorUtil::Add(Velocity, Acceleration);
    
    float DeltaPosition[3];
    VectorUtil::ScalarMultiply(DeltaPosition, Velocity, TimeStep);

    VectorUtil::Add(Position, DeltaPosition);

    // std::cout << PrevPosition[0] << ", " << PrevPosition[1] << ", " << PrevPosition[2] << " --> " << Position[0] << ", " << Position[1] << ", " << Position[2] << std::endl;
}



Simulator::Simulator()
{
	GravitationalParameter = 398600.4418;
    TimeStep = 1.f; // 1.f / 60.f;
	Duration = 10 * 1;
	FieldOfView = 60.f;
	DetectionRange = 1;

	CubesatPosition[0] = 2000.f;
	CubesatPosition[1] = 0.f;
	CubesatPosition[2] = 0.f;

	CubesatVelocity[0] = 0.f;
	CubesatVelocity[1] = 7.5f;
	CubesatVelocity[2] = 0.f;

    DetectionCount = 0;

    NumDebris = -1;
    DebrisPosition = nullptr;
    DebrisVelocity = nullptr;
}


void Simulator::Run()
{
    if (DebrisPosition == nullptr || DebrisVelocity == nullptr) 
    {
        std::cout << "Can't Run Simulation: Debris Position / Velocity unset" << std::endl;
        return;
    }

    DetectionCount = 0;
    int Steps = Duration / TimeStep;

    std::cout << "Start Simulation (steps=" << Steps << ")" << std::endl;

    for (int i = 0; i < Steps; i++) 
    {
        UpdateOrbitalBody(CubesatPosition, CubesatVelocity, GravitationalParameter, TimeStep);
        
        for (int d = 0; d < NumDebris; d++)
        {
            // if mod(i, floor(this.num_debris / 10)) == 0
            //    fprintf("Debris %d\n", i)

            UpdateOrbitalBody(&DebrisPosition[d*3], &DebrisVelocity[d*3], GravitationalParameter, TimeStep);


            float Difference[3];
            VectorUtil::Sub(Difference, &DebrisPosition[d*3], CubesatPosition);

            float DistanceSquared = Difference[0] * Difference[0] + Difference[1] * Difference[1] + Difference[2] * Difference[2];
            // std::cout << DistanceSquared << std::endl;
            bool bDetected = DistanceSquared < (DetectionRange * DetectionRange);

            DetectionCount += bDetected;
        }
        std::cout << "Simulator Iteration" << std::endl;
    }
    std::cout << "Simulator Finished" << std::endl;

    free(DebrisPosition);
    free(DebrisVelocity);
}
