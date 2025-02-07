#pragma once
class DebrisDataGenerator
{
private:
	float* DebrisPosition;
	float* DebrisVelocity;
	int NumDebris;

public:
	DebrisDataGenerator();

public:
	void FetchData();
	void CacheData();
	void LoadData();

public:
	void TakeData(float*& OutDebrisPosition, float*& OutDebrisVelocity, int& OutNumDebris);
};

