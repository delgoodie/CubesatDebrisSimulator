#pragma once

#include <cmath>
#include <glm/glm.hpp>


struct Debris {
    float Ecc = 0.f;
    float Inc = 0.f;
    float RA = 0.f;
    float AP = 0.f;
    float MA = 0.f;
};


struct DebrisList {
    Debris* head;
    size_t num;
};


class VectorUtil
{


public:
    static void Copy(float Target[3], float Source[3]);
    static void ScalarMultiply(float Vec[3], float Scalar);
    static void ScalarMultiply(float Target[3], float Source[3], float Scalar);
    static void Add(float Vec[3], float Addition[3]);
    static void Add(float Target[3], float Source[3], float Addition[3]);
    static void Sub(float Vec[3], float Subtraction[3]);
    static void Sub(float Target[3], float Source[3], float Subtraction[3]);
    static float Size(float Vector[3]);
    static float Dot(float A[3], float B[3]);
    static void Cross(float Out[3], float A[3], float B[3]);
    static void Square(float Vector[3]);
};



void COE_from_SV(float Position[3], float Velocity[3], float mu)
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

    float RA;
    const float PI = 3.14159265;




    float eps = 1e-10;
    float r = VectorUtil::Size(Position);
    float v = VectorUtil::Size(Velocity);

    float vr = VectorUtil::Dot(Position, Velocity) / r;

    float H[3];
    VectorUtil::Cross(H, Position, Velocity);
    float h = VectorUtil::Size(H);

    // Equation 4.7:
    float incl = acos(H[2] / h);

    // Equation 4.8:
    float ZAxis[3] = { 0.f, 0.f, 1.f };
    float N[3];
    VectorUtil::Cross(N, ZAxis, H);
    float n = VectorUtil::Size(N);

    // Equation 4.9:
    if (n != 0)
    {
        RA = acos(N[0] / n);
        if (N[2] < 0)
        {
            RA = 2 * PI - RA;
        }
    }
    else
    {
        RA = 0;
    }

    /*
    // Equation 4.10:
    float E = 1 / mu * ((v ^ 2 - mu / r) * R - r * vr * V);
    float e = norm(E);

    // Equation 4.12 (incorporating the case e = 0) :
    if (n != 0) 
    {

        if (e > eps) 
        {
            w = acos(dot(N, E) / n / e);
            if (E(3) < 0) 
            {
                w = 2 * pi - w;
            }
        }
        else 
        {
            w = 0;
        }
    else
    {
        w = 0;
    }


    // Equation 4.13a(incorporating the case e = 0) :
    if (e > eps) 
    {
        TA = acos(dot(E, R) / e / r);
        if (vr < 0)
        {
            TA = 2 * pi - TA;
        }
    }
    else
    {
        cp = cross(N, R);
        if (cp(3) >= 0)
        {
            TA = acos(dot(N, R) / n / r);
        }
        else
        {
            TA = 2 * pi - acos(dot(N, R) / n / r);
        }
    }
    // Equation 4.62 (a < 0 for a hyperbola) :
    a = h ^ 2 / mu / (1 - e ^ 2);
    coe = [h e RA incl w TA a];

    */
}
