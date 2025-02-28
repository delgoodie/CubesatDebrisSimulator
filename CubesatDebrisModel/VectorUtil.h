#pragma once

#include <cmath>
#include <glm/glm.hpp>


struct CoordKep {
    float a = 0.f;
    float e = 0.f;
    float i = 0.f;
    float o = 0.f;
    float w = 0.f;
    float m = 0.f;
};

struct CoordCar
{
    glm::vec3 pos;
    glm::vec3 vel;
};


struct Debris
{
    CoordKep coord;
    float Size;
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

public:
    static float TA_to_MA(float TA, float Ecc);
    static float MA_to_TA(float MA, float Ecc);


    static CoordKep CC_to_CK(const CoordCar& CC);
    static CoordCar CK_to_CC(const CoordKep& CK);
};