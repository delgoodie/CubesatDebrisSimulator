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
