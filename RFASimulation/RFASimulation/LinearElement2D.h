#pragma once
#include <vector>


namespace RFASimulation
{
	class LinearElement2D
	{
		public: 
			
			std::vector<std::vector<int>> StiffnessMatrix;

			std::vector<std::vector<int>> MassMatrix;

			std::vector<int> RightHandSide;
			
			
			LinearElement2D(void);

			~LinearElement2D(void);


	};

}

