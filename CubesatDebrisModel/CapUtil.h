#pragma once

#include <cmath>
#include <glm/glm.hpp>
#include <stdexcept>


typedef glm::highp_dvec3 vec3;

struct CoordKep {
    double a = 0;
    double e = 0;
    double i = 0;
    double o = 0;
    double w = 0;
    double m = 0;

    bool Equals(const CoordKep& Other) 
    {
        return abs(a - Other.a) < .001 && abs(e - Other.e) < .0001 && abs(i - Other.i) < .0001 && abs(o - Other.o) < .0001 && abs(w - Other.w) < .0001 && abs(m - Other.m) < .0001;
    }
};


struct CoordKepF {
    float a = 0;
    float e = 0;
    float i = 0;
    float o = 0;
    float w = 0;
    float m = 0;
};

struct CoordCar
{
    vec3 pos;
    vec3 vel;
};


struct Debris
{
    CoordKep coord;
    double rcs;
};


struct DebrisF
{
    CoordKepF coord;
    float rcs;
};




struct DebrisList {
private:
    std::unique_ptr<Debris[]> head;
public:
    const size_t num;

public:
    DebrisList() : head(nullptr), num(0) {}
    DebrisList(size_t InNum) : head(std::make_unique<Debris[]>(InNum)), num(InNum) { }
    DebrisList(DebrisList&& other) noexcept;

    DebrisList(const DebrisList&) = delete;
    DebrisList& operator=(const DebrisList&) = delete;

    DebrisList& operator=(DebrisList&& other) noexcept
    {
        if (this != &other)
        {
            head = std::move(other.head);
            *const_cast<size_t*>(&num) = other.num;
            *const_cast<size_t*>(&other.num) = 0;
        }
        return *this;
    }

    ~DebrisList();


    Debris& operator[](size_t index) 
    {
        if (index >= num) throw std::out_of_range("Index out of bounds");
        return head.get()[index];
    }

    const Debris& operator[](size_t index) const 
    {
        if (index >= num) throw std::out_of_range("Index out of bounds");
        return head.get()[index];
    }

    Debris* GetRawData() { return head.get(); }

    bool IsValid() { return num > 0 && head.get() != nullptr; }
};


class CapUtil
{
public:
    static double TA_to_MA(double TA, double Ecc);
    static double MA_to_TA(double MA, double Ecc);


    static CoordKep CC_to_CK(const CoordCar& CC);
    static CoordCar CK_to_CC(const CoordKep& CK);


    static void TestMinDistanceAlgorithm();
    // static double MinDistance(const CoordKep& orbitA, const CoordKep& orbitB, std::vector<vec3>& SamplesA, std::vector<vec3>& SamplesB);
    static double MinDistanceBetweenEllipses(const CoordKep& orbitA, const CoordKep& orbitB);

    static double Deg2Rad(const double& AngleDegrees) { return AngleDegrees * 3.1415926535 / 180.; }
    static double Rad2Deg(const double& AngleRadians) { return AngleRadians * 180. / 3.1415926535; }


    static double FRand(double Min, double Max) 
    {
        return  Min + static_cast <double> (rand()) / (static_cast <double> (RAND_MAX / (Max - Min)));
    }
};