#pragma once
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
};

