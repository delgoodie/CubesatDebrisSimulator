#pragma once
#include <string>
#include "CapUtil.h"

class DataManager
{

public:
	DebrisList GetData();

private:
	void WriteDataToBinFile(const char* FileName, DebrisList& debrisList);
	void WriteDataToCSV(const char* FileName, DebrisList& debrisList);

	DebrisList FetchMocatData_csv(const char* FileName);
	DebrisList FetchLeoData_json(const char* FileName);
	DebrisList FetchData_bin(const char* FileName);

};