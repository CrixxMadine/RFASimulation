#include "TriangleMesh.h"

using namespace RFASimulation;


TriangleMesh::TriangleMesh() 
{
    TriangleMesh::m_tmesh = vector<vector<int>>();
    TriangleMesh::m_pmesh = vector<vector<double>>();
    TriangleMesh::m_bedges = vector<vector<double>>();
}

TriangleMesh::TriangleMesh(vector<vector<int>> tmesh, vector<vector<double>> pmesh, vector<vector<double>> bedges)
{
    TriangleMesh::m_tmesh = tmesh;
    TriangleMesh::m_pmesh = pmesh;
    TriangleMesh::m_bedges = bedges;
}


TriangleMesh::~TriangleMesh(void) {};

int TriangleMesh::GetPointNumber(void)
{
    return m_pmesh[0].size();
}

int TriangleMesh::GetTriangleNumber(void)
{
    return m_tmesh[0].size();
}
