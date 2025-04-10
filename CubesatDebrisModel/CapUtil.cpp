#include "CapUtil.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>


DebrisList::DebrisList(DebrisList&& other) noexcept
    : head(std::move(other.head)), num(other.num)
{
    *const_cast<size_t*>(&other.num) = 0;
    // std::cout << "DebrisList moved num = " << num << std::endl;
}

DebrisList::~DebrisList()
{
    if (num > 0)
    {
        std::cout << "~DebrisList num = " << num << std::endl;
    }
}



double CapUtil::TA_to_MA(double TA, double Ecc)
{
    double E = 2 * atan(sqrt((1 - Ecc) / (1 + Ecc)) * tan(TA / 2));
    return E - Ecc * sin(E);
}

double CapUtil::MA_to_TA(double MA, double Ecc)
{
    const double TOLERANCE = 1e-6;

    double E = MA;
    double delta = 1.0;

    // Newton-Raphson iteration
    while (abs(delta) > TOLERANCE) 
    {
        delta = (E - Ecc * sin(E) - MA) / (1 - Ecc * cos(E));
        E -= delta;
    }

    return 2. * double(std::atan2(sqrt(1 + Ecc) * sin(E / 2), sqrt(1 - Ecc) * cos(E / 2)));
}

CoordKep CapUtil::CC_to_CK(const CoordCar& CC)
{
    /*
    This function computes the classical orbital elements(coe)
    from the state vector(R, V) using Algorithm 4.1.

    mu - gravitational parameter(km ^ 3 / s ^ 2)
    R - position vector in the geocentric equatorial frame(km)
    V - velocity vector in the geocentric equatorial frame(km)
    r, v - the magnitudes of R and V
    vr - radial velocity component(km / s)
    H - the angular momentum vector(km ^ 2 / s)
    h - the magnitude of H(km ^ 2 / s)
    incl - inclination of the orbit(rad)
    N - the node line vector(km ^ 2 / s)
    n - the magnitude of N
    cp - cross product of N and R
    RA - right ascension of the ascending node(rad)
    E - eccentricity vector
    e - eccentricity(magnitude of E)
    eps - a small number below which the eccentricity is considered
    to be zero
    w - argument of perigee(rad)
    TA - true anomaly(rad)
    a - semimajor axis(km)
    pi - 3.1415926
    coe - vector of orbital elements[h e RA incl w TA a]
    User M - functions required : None
    */

    CoordKep OC;


    const double MU = 3.986004418e5f; // km3/s2
    const double PI = 3.14159265358979323846264;


    double eps = 1e-10f;
    double r = glm::length(CC.pos);
    double v = glm::length(CC.vel);

    double vr = glm::dot(CC.pos, CC.vel) / r;

    vec3 H = glm::cross(CC.pos, CC.vel);
    double h = glm::length(H);

    // Equation 4.7:
    OC.i = acos(H.z / h);

    // Equation 4.8:
    vec3 ZAxis = { 0., 0., 1. };
    glm::vec N = glm::cross(ZAxis, H);
    double n = glm::length(N);

    // Equation 4.9:
    if (n != 0)
    {
        OC.o = acos(N.x / n);
        if (N.z < 0)
        {
            OC.o = 2 * PI - OC.o;
        }
    }
    else
    {
        OC.o = 0;
    }

    // Equation 4.10:
    vec3 E = 1 / MU * ((v * v - MU / r) * CC.pos - r * vr * CC.vel);
    OC.e = glm::length(E);

    // Equation 4.12 (incorporating the case e = 0) :
    if (n != 0)
    {

        if (OC.e > eps)
        {
            OC.w = acos(dot(N, E) / n / OC.e);
            if (E.z < 0)
            {
                OC.w = 2. * PI - OC.w;
            }
        }
        else
        {
            OC.w = 0;
        }
    }
    else
    {
        OC.w = 0;
    }


    // Equation 4.13a(incorporating the case e = 0)
    double TA;
    if (OC.e > eps)
    {
        TA = acos(glm::dot(E, CC.pos) / OC.e / r);
        if (vr < 0)
        {
            TA = 2 * PI - TA;
        }
    }
    else
    {
        vec3 cp = glm::cross(N, CC.pos);
        if (cp.z >= 0)
        {
            TA = acos(glm::dot(N, CC.pos) / n / r);
        }
        else
        {
            TA = 2 * PI - acos(glm::dot(N, CC.pos) / n / r);
        }
    }
    // Equation 4.62 (a < 0 for a hyperbola) :
    OC.a = h * h / MU / (1 - OC.e * OC.e);

    OC.m = TA_to_MA(TA, OC.e);

    // didn't save
    // h, a

    return OC;
}

CoordCar CapUtil::CK_to_CC(const CoordKep& CK)
{
    const double MU = 3.986004418e5f; // km3/s2
    const double PI = 3.14159265358979323846264;
    const double TOLERANCE = 1e-6; // Convergence tolerance for Kepler's equation

    CoordCar CC;

    if (CK.e > 1) 
    {
        std::cout << "HELP" << std::endl;
    }


    assert(CK.e <= 1);
    if (CK.e > 1) return {};

    double E = CK.m;
    double delta = 1.0;
    while (abs(delta) > TOLERANCE)
    {
        delta = (E - CK.e * sin(E) - CK.m) / (1 - CK.e* cos(E));
        E -= delta;
    }
    double TA = 2. * double( atan2(sqrt(1 + CK.e) * sin(E / 2), sqrt(1 - CK.e) * cos(E / 2)) );
    
    double r = CK.a * (1 - CK.e * CK.e) / (1 + CK.e * cos(TA));

    // Compute perifocal coordinates (PQW)
    double x_P = r * cos(TA);
    double y_P = r * sin(TA);
    double z_P = 0.0;

    // Compute mean motion n
    double n = std::sqrt(MU / (CK.a * CK.a * CK.a));

    // Compute velocity in perifocal coordinates
    double v_x = -sin(E) * n * CK.a / (1 - CK.e * cos(E));
    double v_y = sqrt(1 - CK.e * CK.e) * cos(E) * n * CK.a / (1 - CK.e * cos(E));
    double v_z = 0.0;

    // Rotation matrix from perifocal to inertial coordinates
    double cosO = cos(CK.o), sinO = sin(CK.o);
    double coso = cos(CK.o), sino = sin(CK.o);
    double cosi = cos(CK.i), sini = sin(CK.i);

    // Position transformation
    CC.pos.x = double( (cosO * coso - sinO * sino * cosi) * x_P + (-cosO * sino - sinO * coso * cosi) * y_P );
    CC.pos.y = double( (sinO * coso + cosO * sino * cosi) * x_P + (-sinO * sino + cosO * coso * cosi) * y_P );
    CC.pos.z = double( (sino * sini) * x_P + (coso * sini) * y_P );

    // Velocity transformation
    CC.vel.x = double( (cosO * coso - sinO * sino * cosi) * v_x + (-cosO * sino - sinO * coso * cosi) * v_y );
    CC.vel.y = double( (sinO * coso + cosO * sino * cosi) * v_x + (-sinO * sino + cosO * coso * cosi) * v_y );
    CC.vel.z = double( (sino * sini) * v_x + (coso * sini) * v_y );

    return CC;
}



void CapUtil::TestMinDistanceAlgorithm()
{

    double OverallMin = 0;
    CoordKep BestCoordA, BestCoordB;

    std::vector<vec3> SamplesA;
    std::vector<vec3> SamplesB;

    for (int i = 0; i < 10000; i++)
    {
        SamplesA.clear();
        SamplesB.clear();

        CoordKep CoordA = { 8000., 0., CapUtil::Deg2Rad(CapUtil::FRand(0., 90.)), CapUtil::Deg2Rad(CapUtil::FRand(0., 90.)), CapUtil::Deg2Rad(CapUtil::FRand(0., 90.)), 0. };
        CoordKep CoordB = { 8000., 0., CapUtil::Deg2Rad(CapUtil::FRand(0., 90.)), CapUtil::Deg2Rad(CapUtil::FRand(0., 90.)), CapUtil::Deg2Rad(CapUtil::FRand(0., 90.)), 0. };

        double MinDist = CapUtil::MinDistanceBetweenEllipses(
            CoordA, CoordB
            // SamplesA,
            // SamplesB
        );
        if (MinDist > OverallMin) {
            OverallMin = MinDist;
            BestCoordA = CoordA;
            BestCoordB = CoordB;
        }
    }


    std::cout << OverallMin << std::endl;
    std::cout << CapUtil::Rad2Deg(BestCoordA.i) << "  " << CapUtil::Rad2Deg(BestCoordA.o) << "  " << CapUtil::Rad2Deg(BestCoordA.w) << std::endl;
    std::cout << CapUtil::Rad2Deg(BestCoordB.i) << "  " << CapUtil::Rad2Deg(BestCoordB.o) << "  " << CapUtil::Rad2Deg(BestCoordB.w) << std::endl;
    SamplesA.clear();
    SamplesB.clear();
    double MinDist2 = CapUtil::MinDistanceBetweenEllipses(
        BestCoordA, BestCoordB
        // SamplesA,
        // SamplesB
    );


    std::ofstream outFile("out_test.csv");

    if (!outFile) {
        std::cerr << "Error opening file for writing!\n";
    }

    outFile << "Ax, Ay, Az, Bx, By, Bz\n";

    for (int i = 0; i < SamplesA.size(); i++)
    {
        outFile << SamplesA[i].x << ", " << SamplesA[i].y << ", " << SamplesA[i].z << ", " << SamplesB[i].x << ", " << SamplesB[i].y << ", " << SamplesB[i].z << "\n";
    }

    outFile.close();
    std::cout << "CSV file written successfully!\n";

}

// double CapUtil::MinDistance(const CoordKep& orbitA, const CoordKep& orbitB, std::vector<vec3>& SamplesA, std::vector<vec3>& SamplesB)
double CapUtil::MinDistanceBetweenEllipses(const CoordKep& orbitA, const CoordKep& orbitB)
{
    const double PI = 3.14159265358979323846264;

    double best_mA = 0., best_mB = 0.;
    double minDistSqr = std::numeric_limits<double>::max();

    std::vector<double> IterWindow = { 2. * PI, PI / 6., PI / 60., PI / 200., PI / 800., PI / 1000., PI / 2000., PI / 15000., PI / 60000., PI / 300000. };

    for (int Iter = 0; Iter < IterWindow.size(); Iter++)
    {
        double min_mA = best_mA - IterWindow[Iter] / 2.;
        double max_mA = best_mA + IterWindow[Iter] / 2.;
        double min_mB = best_mB - IterWindow[Iter] / 2.;
        double max_mB = best_mB + IterWindow[Iter] / 2.;


        const int samples = 16;
        minDistSqr = std::numeric_limits<double>::max();

        std::vector<std::tuple<vec3, double>> pointsA;
        for (int i = 0; i <= samples; ++i)
        {
            CoordKep sampleA = orbitA;
            sampleA.m = min_mA + (max_mA - min_mA) * (double)i / (double)samples;
            vec3 posA = CapUtil::CK_to_CC(sampleA).pos;
            pointsA.push_back({ posA, sampleA.m });
            // SamplesA.push_back(posA);
        }
        std::vector<std::tuple<vec3, double>> pointsB;
        for (int i = 0; i <= samples; ++i)
        {
            CoordKep sampleB = orbitB;
            sampleB.m = min_mB + (max_mB - min_mB) * (double)i / (double)samples;
            vec3 posB = CapUtil::CK_to_CC(sampleB).pos;
            pointsB.push_back({ posB, sampleB.m });
            // SamplesB.push_back(posB);
        }

        for (const std::tuple<vec3, double>& posA : pointsA)
        {
            for (const std::tuple<vec3, double>& posB : pointsB)
            {
                vec3 dif = std::get<vec3>(posA) - std::get<vec3>(posB);
                double distSqr = glm::dot(dif, dif);

                if (distSqr < minDistSqr)
                {
                    minDistSqr = distSqr;
                    best_mA = std::get<double>(posA);
                    best_mB = std::get<double>(posB);
                }
            }
        }

        
    }
    return sqrt(minDistSqr);
}