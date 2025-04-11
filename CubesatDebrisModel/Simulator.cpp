#include "Simulator.h"

#include <iostream>
#include "CapUtil.h"
#include <vector>
#include <chrono>
#include <random>


using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;


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

    auto total_t1 = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    std::vector<Debris> DebrisInRange;
    for (int i = 0; i < debrisList.num; i++) 
    {
        if ((i % int(debrisList.num / 100)) == 0 && i > 0)
        {
            auto t2 = high_resolution_clock::now();
            auto ms_int = duration_cast<milliseconds>(t2 - t1);
            duration<double, std::milli> ms_double = t2 - t1;

            int Percent = i / (debrisList.num / 100.);

            std::cout << "Cull Iter debris: " << i << "  " << Percent << "%   " << double(ms_int.count())/1000. << " s " << std::endl;
            t1 = high_resolution_clock::now();
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

    auto total_t2 = high_resolution_clock::now();
    auto total_ms_int = duration_cast<milliseconds>(total_t2 - total_t1);
    std::cout << "Culling Finished " << InitialSize - FinalSize << " debris out of " << InitialSize << " initial debris (" << (1. - ((double)FinalSize / double(InitialSize))) * 100. << "%)   in " << double(total_ms_int.count()) / 1000. << "s" << std::endl;
}

void Simulator::Run()
{
    if (debrisList.num == 0) 
    {
        std::cout << "Can't Run Simulation: Debris Position / Velocity unset" << std::endl;
        return;
    }

    DetectionCount = 0;
    
    auto total_t1 = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    for (int i = 0; i < debrisList.num; i++) 
    {
        if ((i % (debrisList.num / 100)) == 0 && i > 0)
        {
            auto t2 = high_resolution_clock::now();
            auto ms_int = duration_cast<milliseconds>(t2 - t1);
            duration<double, std::milli> ms_double = t2 - t1;

            int Percent = i / (debrisList.num / 100.);
            
            std::cout << "Sim Iter debris: " << i << "  " << Percent << "%   " << ms_int.count() << " ms " << std::endl;
            t1 = high_resolution_clock::now();
        }

        SimulateDebrisCubesatPair(debrisList[i]);
    }

    auto total_t2 = high_resolution_clock::now();
    auto total_ms_int = duration_cast<milliseconds>(total_t2 - total_t1);
    std::cout << "Sim Finished  detected: " << DetectionCount << "  in " << double(total_ms_int.count()) / 1000. << "s" << std::endl;
}

void Simulator::SimulateDebrisCubesatPair(const Debris& debris)
{
    CoordKep db = debris.coord;
    CoordKep cs = cubesat.coord;
    double mOffset = db.m;

    const double noise_mean = .1;
    const double noise_stdev = .05;
    const double signal_std_min = 2;
    const double acc_signal_threshold = .01;
    double signal_acc = 0.0;

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(noise_mean, noise_stdev);


    for (double acc_m = 0.; acc_m < 2 * PI;)
    {
        double prev_m = cs.m;

        CoordCar dcc = CapUtil::CK_to_CC(db);
        CoordCar csc = CapUtil::CK_to_CC(cs);
        double rcs = 10.;

        double distance = glm::distance(dcc.pos, csc.pos);

        if (distance < cubesat.DetectionRange * 2)
        {
            const double TimeStep = .001;

            double noise = distribution(generator);
            double echo = pow(cubesat.TransmitPower * cubesat.TransmitAntennaGainFactor * cubesat.TransmitAntennaGainFactor * cubesat.Wavelength * cubesat.Wavelength * rcs / pow(4*PI, 3), .25);

            double signal = echo + noise;

            double strength = signal - (noise_mean + noise_stdev * signal_std_min);
            if (strength > 0.) 
            {
                signal_acc += strength * TimeStep;
            }

            if (signal_acc >= acc_signal_threshold)
            {
                DetectionCount++;
                return;
            }

            UpdateOrbitalBody(db, GravitationalParameter, TimeStep);
            UpdateOrbitalBody(cs, GravitationalParameter, TimeStep);
        }
        else
        {
            vec3 normDir = (dcc.pos - csc.pos) / distance;
            double COM_vel = glm::dot(csc.vel, normDir) - glm::dot(dcc.vel, normDir);

            double TimeStep = 500.;

            if (COM_vel > 1.0)
            {
                double minTime = distance / COM_vel;

                if (minTime < 500.0)
                {
                    TimeStep = minTime / 3.;
                }
            }

            UpdateOrbitalBody(db, GravitationalParameter, TimeStep);
            UpdateOrbitalBody(cs, GravitationalParameter, TimeStep);
        }

        double delta_m = cs.m - prev_m;
        if (delta_m > PI) delta_m = 2. * PI - delta_m;
        acc_m += abs(delta_m);
    }
}


