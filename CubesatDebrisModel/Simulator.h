#pragma once




class Simulator
{
private:
	float GravitationalParameter; // (mu) km3/s2
	float TimeStep; // s
	float Duration;
	float FieldOfView;
	float DetectionRange;

	int NumDebris;

	float* DebrisPosition;
	float* DebrisVelocity;

	float CubesatPosition[3];
	float CubesatVelocity[3];

public:
	int DetectionCount;


public:
	Simulator();

public:
	void SetDebris(float* InDebrisPosition, float* InDebrisVelocity, int InNumDebris) { DebrisPosition = InDebrisPosition; DebrisVelocity = InDebrisVelocity; NumDebris = InNumDebris; }
	void Run();

};

