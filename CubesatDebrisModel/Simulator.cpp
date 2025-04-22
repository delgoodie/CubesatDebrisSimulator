#include "Simulator.h"

#include <iostream>
#include "CapUtil.h"
#include <vector>
#include <chrono>
#include <random>
#include <deque>
#include <algorithm>


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
        if (debrisList.num > 100 && (i % int(debrisList.num / 100)) == 0 && i > 0)
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
void Simulator::CullDebrisByIntersection()
{
    std::cout << "Culling debris by intersection" << std::endl;

    int InitialSize = (int)debrisList.num;

    auto total_t1 = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    std::vector<Debris> DebrisInRange;
    for (int i = 0; i < debrisList.num; i++)
    {
        if (debrisList.num > 100 && (i % int(debrisList.num / 100)) == 0 && i > 0)
        {
            auto t2 = high_resolution_clock::now();
            auto ms_int = duration_cast<milliseconds>(t2 - t1);
            duration<double, std::milli> ms_double = t2 - t1;

            int Percent = i / (debrisList.num / 100.);

            std::cout << "Cull Iter debris: " << i << "  " << Percent << "%   " << double(ms_int.count()) / 1000. << " s " << std::endl;
            t1 = high_resolution_clock::now();
        }

        if (FindDebrisCubesatPairIntersectionTime(debrisList[i]).has_value())
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

    int NumInRange = 0;
    DetectionCount = 0;
    
    auto total_t1 = high_resolution_clock::now();
    auto t1 = high_resolution_clock::now();

    for (int i = 0; i < debrisList.num; i++) 
    {
        if (debrisList.num > 100 && (i % (debrisList.num / 100)) == 0 && i > 0)
        {
            auto t2 = high_resolution_clock::now();
            auto ms_int = duration_cast<milliseconds>(t2 - t1);
            duration<double, std::milli> ms_double = t2 - t1;

            int Percent = i / (debrisList.num / 100.);
            
            std::cout << "Sim Iter debris: " << i << "  " << Percent << "%   " << ms_int.count() << " ms " << std::endl;
            t1 = high_resolution_clock::now();
        }

        std::optional<double> IntersectionStartTime = FindDebrisCubesatPairIntersectionTime(debrisList[i]);

        if (IntersectionStartTime.has_value()) 
        {
            NumInRange++;
            SimulateDebrisCubesatPair(debrisList[i], IntersectionStartTime.value());
        }
    }

    int MissedInRange = NumInRange - DetectionCount;

    auto total_t2 = high_resolution_clock::now();
    auto total_ms_int = duration_cast<milliseconds>(total_t2 - total_t1);
    std::cout << "Sim Finished  detected: " << DetectionCount << "  in " << double(total_ms_int.count()) / 1000. << "s      missed = " << MissedInRange << std::endl;
}

std::optional<double> Simulator::FindDebrisCubesatPairIntersectionTime(const Debris& debris)
{
    CoordKep db = debris.coord;
    CoordKep cs = cubesat.coord;
    double mOffset = db.m;
    double MaxTimeStep = 500.0;
    const double MinTimeStep = .0001;
    double StabilityFactor = 3.0;

    const double Tolerance = cubesat.DetectionRange / 100.0; // get within 1% of start of detection range

    double Time = 0.0;
    for (double acc_m = 0.; acc_m < 2 * PI;)
    {
        double prev_m = cs.m;

        CoordCar dcc = CapUtil::CK_to_CC(db);
        CoordCar csc = CapUtil::CK_to_CC(cs);

        double distance = glm::distance(dcc.pos, csc.pos);

        vec3 normDir = (dcc.pos - csc.pos) / distance;
        double COM_vel = glm::dot(csc.vel, normDir) - glm::dot(dcc.vel, normDir);
        double displacement = distance - cubesat.DetectionRange;

        if (displacement < 0.0 && abs(displacement / COM_vel) < MinTimeStep)
        {
            std::cout << displacement * 1000 << " m displacement" << std::endl;
            return Time;
        }

        double TimeStep = MaxTimeStep;
           
        if (COM_vel > 1.0)
        {
            double MinTimeToIntersection = displacement / COM_vel;
            TimeStep = std::clamp(MinTimeToIntersection / StabilityFactor, MinTimeStep, MaxTimeStep);
        }

        UpdateOrbitalBody(db, GravitationalParameter, TimeStep);
        UpdateOrbitalBody(cs, GravitationalParameter, TimeStep);
         
        Time += TimeStep;
       
        double delta_m = cs.m - prev_m;
        if (delta_m > PI) delta_m = 2. * PI - delta_m;
        acc_m += abs(delta_m);
    }

    return std::optional<double>();
}

void Simulator::SimulateDebrisCubesatPair(const Debris& debris, double StartTime)
{
    const double BoltzmansConstant = 1.38e-23;  // Boltzmann constant (J/K)
    const double SystemTemperature = 10.0;  // System noise temperature (K)
    const double Bandwidth = 1e6;       // Bandwidth (Hz)
    const double TimeStep = .0001;
    const double NoiseMean = BoltzmansConstant * SystemTemperature * Bandwidth;
    const double RadarCrossSection = CapUtil::FRand(0.0001, 0.00001); // m2
    const int WindowSize = 100;
    const double DetectionThresholdPower = WindowSize * NoiseMean + 4.0 * NoiseMean * sqrt(WindowSize); // noise sum mean plus 4 sigma

    std::default_random_engine generator;
    std::exponential_distribution<double> NoiseDistribution(1.0 / NoiseMean);

    CoordKep db = debris.coord;
    CoordKep cs = cubesat.coord;

    // fast forward to Start Time
    UpdateOrbitalBody(db, GravitationalParameter, StartTime);
    UpdateOrbitalBody(cs, GravitationalParameter, StartTime);

    CoordCar dcc = CapUtil::CK_to_CC(db);
    CoordCar csc = CapUtil::CK_to_CC(cs);
    
    std::deque<double> PowerReceivedWindow;
    double IntegratedPower = 0.0;
    for (int i = 0; i < WindowSize; i++)
    {
        double PowerNoise = NoiseDistribution(generator);
        IntegratedPower += PowerNoise;
        PowerReceivedWindow.push_back(PowerNoise);
    }

    double Distance = glm::distance(dcc.pos, csc.pos);

    while (Distance < cubesat.DetectionRange)
    {
        double PowerNoise = NoiseDistribution(generator);

        double PowerEcho = (cubesat.TransmitPower * pow(cubesat.TransmitAntennaGainFactor, 2) * pow(cubesat.Wavelength, 2) * RadarCrossSection)
            /
            (pow(4 * PI, 3) * pow(Distance, 4) * cubesat.SystemLoss);

        double PowerReceived = PowerEcho + PowerNoise;

        IntegratedPower += PowerReceived;
        IntegratedPower -= PowerReceivedWindow[0];

        PowerReceivedWindow.push_back(PowerReceived);
        PowerReceivedWindow.pop_front();

        if (IntegratedPower > DetectionThresholdPower)
        {
            DetectionCount++;
            return;
        }

        // Update for next frame

        UpdateOrbitalBody(db, GravitationalParameter, TimeStep);
        UpdateOrbitalBody(cs, GravitationalParameter, TimeStep);

        dcc = CapUtil::CK_to_CC(db);
        csc = CapUtil::CK_to_CC(cs);

        Distance = glm::distance(dcc.pos, csc.pos);
    }
}



double computeSunVisibility(const vec3& r_sat, const vec3& u_sun) 
{
    double R_e = 6378.0;
    double h = glm::dot(r_sat, u_sun);

    if (h > 0) 
    {
        return 1.0;
    }

    double r2 = dot(r_sat, r_sat);
    double d2 = r2 - h * h;

    if (d2 < R_e * R_e) 
    {
        return 0.0;
    }
    else 
    {
        return 1.0;
    }
}




std::vector<double> Simulator::GenerateTemperatureField(double TimeStep)
{
    int Day = 1;

    //--- constants & parameters
    const double Period = 2 * PI * sqrt(pow(cubesat.coord.a, 3) / GravitationalParameter);
    const double mass = 1.33; // kg
    const double cp = 900; // J/(kg K)
    const double C = mass * cp; // J/K
    const double alpha = 0.7; // solar absorptivity
    const double eps = 0.8; // IR emissivity
    const double I_sun = 1361.0; // W/m2
    const double rho_alb = 0.3; // Earth albedo
    const double E_IR = 237.0; // W/m2
    const double sigma = 5.670374e-8; // W/m2 K4
    const double side = 0.1; // 10 cm cube
    const double A_total = 6 * side * side; // m2
    const double A_abs = A_total / 4.0; // avg proj area m2
    const double P_int = 2.0; // W, internal electronics

    std::vector<double> Temperatures(int(Period / TimeStep));

    CoordKep cs = cubesat.coord;
    Temperatures.push_back(293.0); // initial temp (K)

    for (double Time = 0.0; Time <= Period; Time += TimeStep)
    {
        UpdateOrbitalBody(cs, GravitationalParameter, TimeStep);

        CoordCar csc = CapUtil::CK_to_CC(cubesat.coord);

        double OldTemp = Temperatures[Temperatures.size() - 1];

        double V_sun = computeSunVisibility(orbit, TimeStep);
        double V_alb = computeAlbVisibility(orbit, TimeStep);
        double V_IR = computeIRVisibility(orbit, TimeStep);

        // heat inputs
        double Q_sun = alpha * A_abs * I_sun * V_sun;
        double Q_alb = alpha * A_abs * I_sun * rho_alb * V_alb;
        double Q_ir = eps * A_abs * E_IR * V_IR;
        double Q_int = P_int;

        // radiative loss
        double Q_out = eps * sigma * A_total * pow(OldTemp, 4);

        // temperature update
        double dT = (TimeStep / C) * (Q_sun + Q_alb + Q_ir + Q_int - Q_out);
        double NewTemperature = OldTemp + dT;
        Temperatures.push_back(NewTemperature);
    }

    return Temperatures;
}


