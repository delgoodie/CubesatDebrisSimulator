#include "VectorUtil.h"


void VectorUtil::Copy(float Target[3], float Source[3])
{
    Target[0] = Source[0];
    Target[1] = Source[1];
    Target[2] = Source[2];
}

void VectorUtil::ScalarMultiply(float Vec[3], float Scalar)
{
    Vec[0] *= Scalar;
    Vec[1] *= Scalar;
    Vec[2] *= Scalar;
}

void VectorUtil::ScalarMultiply(float Target[3], float Source[3], float Scalar)
{
    Target[0] = Source[0] * Scalar;
    Target[1] = Source[1] * Scalar;
    Target[2] = Source[2] * Scalar;
}

void VectorUtil::Add(float Vec[3], float Addition[3])
{
    Vec[0] += Addition[0];
    Vec[1] += Addition[1];
    Vec[2] += Addition[2];
}

void VectorUtil::Add(float Target[3], float Source[3], float Addition[3])
{
    Target[0] = Source[0] + Addition[0];
    Target[1] = Source[1] + Addition[1];
    Target[2] = Source[2] + Addition[2];
}

void VectorUtil::Sub(float Vec[3], float Subtraction[3])
{
    Vec[0] -= Subtraction[0];
    Vec[1] -= Subtraction[1];
    Vec[2] -= Subtraction[2];
}

void VectorUtil::Sub(float Target[3], float Source[3], float Subtraction[3])
{
    Target[0] = Source[0] - Subtraction[0];
    Target[1] = Source[1] - Subtraction[1];
    Target[2] = Source[2] - Subtraction[2];
}

float VectorUtil::Size(float Vector[3])
{
    return sqrt(Vector[0] * Vector[0] + Vector[1] * Vector[1] + Vector[2] * Vector[2]);
}

float VectorUtil::Dot(float A[3], float B[3])
{
    return A[0]*B[0] + A[1]*B[1] + A[2]*B[2];
}

void VectorUtil::Cross(float Out[3], float A[3], float B[3])
{
    Out[0] = A[1] * B[2] - A[2] * B[1];
    Out[1] = A[2] * B[0] - A[0] * B[2];
    Out[2] = A[0] * B[1] - A[1] * B[0];
}

void VectorUtil::Square(float Vector[3])
{
    Vector[0] *= Vector[0];
    Vector[1] *= Vector[1];
    Vector[2] *= Vector[2];
}


float VectorUtil::TA_to_MA(float TA, float Ecc)
{
    float E = 2 * atan(sqrt((1 - Ecc) / (1 + Ecc)) * tan(TA / 2);
    return E - Ecc * sin(E);
}

float VectorUtil::MA_to_TA(float MA, float Ecc)
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

    return 2.f * std::atan2(sqrt(1 + Ecc) * sin(E / 2), sqrt(1 - Ecc) * cos(E / 2));
}

CoordKep VectorUtil::CC_to_CK(const CoordCar& CC)
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


    const float MU = 3.986004418e5; // km3/s2
    const float PI = 3.14159265358979323846264;


    float eps = 1e-10;
    float r = CC.pos.length();
    float v = CC.vel.length();

    float vr = glm::dot(CC.pos, CC.vel) / r;

    glm::vec3 H = glm::cross(CC.pos, CC.vel);
    float h = H.length();

    // Equation 4.7:
    OC.i = acos(H.z / h);

    // Equation 4.8:
    glm::vec3 ZAxis = { 0.f, 0.f, 1.f };
    glm::vec N = glm::cross(ZAxis, H);
    float n = N.length();

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
    glm::vec3 E = 1 / MU * ((v * v - MU / r) * CC.pos - r * vr * CC.vel);
    OC.e = E.length();

    // Equation 4.12 (incorporating the case e = 0) :
    if (n != 0)
    {

        if (OC.e > eps)
        {
            OC.w = acos(dot(N, E) / n / OC.e);
            if (E.z < 0)
            {
                OC.w = 2.f * PI - OC.w;
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
    float TA;
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
        glm::vec3 cp = glm::cross(N, CC.pos);
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

CoordCar VectorUtil::CK_to_CC(const CoordKep& CK)
{
    const float MU = 3.986004418e5; // km3/s2
    const float PI = 3.14159265358979323846264;
    const double TOLERANCE = 1e-6; // Convergence tolerance for Kepler's equation

    CoordCar CC;

    double E = CK.m;
    double delta = 1.0;
    while (abs(delta) > TOLERANCE)
    {
        delta = (E - CK.e * sin(E) - CK.m) / (1 - CK.e* cos(E));
        E -= delta;
    }
    float TA = 2.f * atan2(sqrt(1 + CK.e) * sin(E / 2), sqrt(1 - CK.e) * cos(E / 2));
    
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
    CC.pos.x = (cosO * coso - sinO * sino * cosi) * x_P + (-cosO * sino - sinO * coso * cosi) * y_P;
    CC.pos.y = (sinO * coso + cosO * sino * cosi) * x_P + (-sinO * sino + cosO * coso * cosi) * y_P;
    CC.pos.z = (sino * sini) * x_P + (coso * sini) * y_P;

    // Velocity transformation
    CC.vel.x = (cosO * coso - sinO * sino * cosi) * v_x + (-cosO * sino - sinO * coso * cosi) * v_y;
    CC.vel.y = (sinO * coso + cosO * sino * cosi) * v_x + (-sinO * sino + cosO * coso * cosi) * v_y;
    CC.vel.z = (sino * sini) * v_x + (coso * sini) * v_y;

    return CC;
}
