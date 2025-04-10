#include "Simulator.h"

#include <iostream>
#include "CapUtil.h"
#include <vector>

#define PI 3.1415926535

void UpdateOrbitalBody(CoordCar& coord, double GravitationalParameter, double TimeStep)
{
    vec3 PrevPosition = coord.pos;
    
    double Radius = glm::length(coord.pos); // sqrt(Position[0] * Position[0] + Position[1] * Position[1] + Position[2] * Position[2]);

    vec3 Acceleration = coord.pos * -GravitationalParameter / (Radius * Radius * Radius);
    coord.vel += Acceleration * TimeStep;
    coord.pos += coord.vel *TimeStep;

    // std::cout << PrevPosition[0] << ", " << PrevPosition[1] << ", " << PrevPosition[2] << " --> " << Position[0] << ", " << Position[1] << ", " << Position[2] << std::endl;
}


void UpdateOrbitalBody(CoordKep& coord, double GravitationalParameter, double TimeStep)
{
    coord.m = fmod(coord.m + sqrt(GravitationalParameter / (coord.a * coord.a * coord.a)) * TimeStep, 2 * PI);
}



Simulator::Simulator()
{
	GravitationalParameter = 398600.4418;
    TimeStep = .1; // 1. / 60.;
	Duration = 100 * 1;

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
        double MinDist = CapUtil::MinDistanceBetweenEllipses(cubesat.coord, debrisList[i].coord);
        if (MinDist <= cubesat.DetectionRange * 1.5)
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

    std::cout << "Culled " << InitialSize - FinalSize << " debris out of " << InitialSize << " initial debris (" << (1. - ((double)FinalSize / double(InitialSize))) * 100. << "%)" << std::endl;
}

void Simulator::Run()
{
    // CullDebrisByMinDistance();


    if (debrisList.num == 0) 
    {
        std::cout << "Can't Run Simulation: Debris Position / Velocity unset" << std::endl;
        return;
    }

    // DetectionCount = 0;
    // int Steps = int(Duration / TimeStep);
    // std::cout << "Start Simulation (steps=" << Steps << ")" << std::endl;
    
    for (int i = 0; i < debrisList.num; i++) 
    {
        if ((i % (debrisList.num / 100)) == 0)
        {
            std::cout << "debris " << i << std::endl;
        }

        CoordKep dbInitial = debrisList[i].coord;
        CoordKep db = debrisList[i].coord;
        CoordKep cs = cubesat.coord;
        double mOffset = db.m;

        bool bDetected = false;
        for (double acc_m = 0.; acc_m < 2 * PI;)
        {
            double prev_m = cs.m;

            CoordCar dcc = CapUtil::CK_to_CC(db);
            CoordCar csc = CapUtil::CK_to_CC(cs);

            // CoordKep dck = CapUtil::CC_to_CK(dbc);
            // CoordKep csk = CapUtil::CC_to_CK(csc);

            
            double OuterRange = 100.0;
            double OuterTimeStep = 5.;
            double InnerTimeStep = .001;

            double distance = glm::distance(dcc.pos, csc.pos);
            if (distance > OuterRange)
            {
                UpdateOrbitalBody(db, GravitationalParameter, OuterTimeStep);
                UpdateOrbitalBody(cs, GravitationalParameter, OuterTimeStep);
            }
            else if (distance > cubesat.DetectionRange * 2)
            {
                double alpha = (distance - cubesat.DetectionRange * 2.) / (OuterRange - cubesat.DetectionRange * 2.);
                double timestep = (OuterTimeStep / 2.0 - InnerTimeStep) * alpha + InnerTimeStep;
                UpdateOrbitalBody(db, GravitationalParameter, timestep);
                UpdateOrbitalBody(cs, GravitationalParameter, timestep);
            }
            else
            {
                if (distance > cubesat.DetectionRange)
                {
                    std::cout << "DETECTED --- " << i << std::endl;
                    DetectionCount++;
                    bDetected = true;
                    break;
                }

                UpdateOrbitalBody(db, GravitationalParameter, InnerTimeStep);
                UpdateOrbitalBody(cs, GravitationalParameter, InnerTimeStep);
            }

            double delta_m = cs.m - prev_m;
            if (delta_m > PI) delta_m = 2. * PI - delta_m;
            acc_m += abs(delta_m);
        }
    }

    /*
    for (int i = 0; i < Steps; i++) 
    {
        if (i % (Steps / 100) == 0)
        {
            std::cout << "Iter " << i << std::endl;
        }

        UpdateOrbitalBody(cubesat.coord, GravitationalParameter, TimeStep);
        vec3 cubesatPos = CapUtil::CK_to_CC(cubesat.coord).pos;

        
        for (int d = 0; d < debrisList.num; d++)
        {
            // if mod(i, floor(this.num_debris / 10)) == 0
            //    fprintf("Debris %d\n", i)

            UpdateOrbitalBody(debrisList[d].coord, GravitationalParameter, TimeStep);

            vec3 debrisPos = CapUtil::CK_to_CC(debrisList[i].coord).pos;

            vec3 difference = cubesatPos - debrisPos;

            double DistanceSquared = glm::dot(difference, difference);
            // std::cout << DistanceSquared << std::endl;
            bool bDetected = DistanceSquared < (cubesat.DetectionRange * cubesat.DetectionRange);

            DetectionCount += bDetected;
        }
    }
    */
    std::cout << "Simulator Finished" << std::endl;
}
