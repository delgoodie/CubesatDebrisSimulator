#pragma once

#include <cmath>
#include <glm/glm.hpp>
#include <stdexcept>

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
    static float TA_to_MA(float TA, float Ecc);
    static float MA_to_TA(float MA, float Ecc);


    static CoordKep CC_to_CK(const CoordCar& CC);
    static CoordCar CK_to_CC(const CoordKep& CK);

    static float MinDistance(const CoordKep& orbitA, const CoordKep& orbitB, std::vector<glm::vec3>& SamplesA, std::vector<glm::vec3>& SamplesB);

    static float Deg2Rad(const float& AngleDegrees) { return AngleDegrees * 3.1415926535f / 180.f; }
    static float Rad2Deg(const float& AngleRadians) { return AngleRadians * 180.f / 3.1415926535f; }
};