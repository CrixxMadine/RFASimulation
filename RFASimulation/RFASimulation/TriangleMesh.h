#pragma once
#include <vector>

using namespace std;

namespace RFASimulation
{
	class TriangleMesh
	{
	public:

		vector<vector<double>> pmesh;

		vector<vector<int>> tmesh;

		vector<int> bedges;


		TriangleMesh(void);

		~TriangleMesh(void);


	};

}
