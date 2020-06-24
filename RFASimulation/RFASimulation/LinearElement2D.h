#pragma once
#include <vector>

using namespace std;

namespace RFASimulation
{
	class LinearElement2D
	{
		public: 
			
			vector<vector<int>> StiffnessMatrix;

			vector<vector<int>> MassMatrix;

			vector<int> RightHandSide;
			
			
			LinearElement2D(void);

			~LinearElement2D(void);


	};

}

