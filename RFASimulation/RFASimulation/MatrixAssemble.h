#pragma once
#include <vector>

using namespace std;

namespace RFASimulation
{
	static class MatrixAssemble
	{
		public:

			vector<vector<int>> AssembleStiffnessMatrix;

			vector<vector<int>> AssembleMassMatrix;

			vector<vector<int>> AssembleRightHandSide;
	};
}
