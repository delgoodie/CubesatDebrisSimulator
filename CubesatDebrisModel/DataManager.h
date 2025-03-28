#pragma once
#include <string>
#include "CapUtil.h"

class DataManager
{

public:
	void WriteData_bin(const char* FileName, DebrisList& debrisList);
	void WriteData_csv(const char* FileName, DebrisList& debrisList);
	DebrisList FetchMocatData_csv(const char* FileName);
	DebrisList FetchLeoData_json(const char* FileName);
	DebrisList FetchData_bin(const char* FileName);
	DebrisList GenData_Leo(int Num);
	bool HasDataFile(const char* FileName);
};