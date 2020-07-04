#include "TriangleMesh.h"

using namespace RFASimulation;


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
