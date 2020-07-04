#pragma once
#include <vector>

using namespace std;

namespace RFASimulation
{
	static class MatrixAssemble
	{
		public:

			vector<vector<int>> AssembleStiffnessMatrix(void);

			vector<vector<int>> AssembleMassMatrix(void);

			vector<vector<int>> AssembleRightHandSide(void);
	};
}
