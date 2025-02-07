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
