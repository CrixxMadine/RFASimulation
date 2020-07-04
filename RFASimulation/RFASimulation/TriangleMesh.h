#pragma once
#include <vector>

using namespace std;

namespace RFASimulation
{
	class TriangleMesh
	{

	public:

		vector<vector<int>> m_tmesh;

		vector<vector<double>> m_pmesh;

		vector<vector<double>> m_bedges;


		TriangleMesh(vector<vector<int>> tmesh, vector<vector<double>> pmesh, vector<vector<double>> bedges);

		~TriangleMesh(void);


		int GetPointNumber(void);

		int GetTriangleNumber(void);

	};

}
