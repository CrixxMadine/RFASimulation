#pragma once
#include "LinearElement2D.h"

namespace RFASimulation
{
	static class ElementBuilder
	{
		public:

			LinearElement2D GetElement(void);

		private:

			vector<vector<int>,vector<int>> GetStiffnessMatrix(void);

			void GetMassMatrix(void);
	};

}

