#include "../include/CellSurvival.hpp"

//----------------------------
CellSurvival::CellSurvival(std::string cellLine_in, std::string cellGeometryType_in)
{
    cellLine = cellLine_in;
    cellGeometryType = cellGeometryType_in;
}

//----------------------------
void CellSurvival::AddCellSurvivalData(std::vector<std::tuple<double,double,double>>  dataCellSurvival_in)
{
    dataCellSurvival = dataCellSurvival_in;
}