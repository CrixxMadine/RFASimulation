#pragma once
#include "LinearElement2D.h"

namespace RFASimulation
{
	static class ElementBuilder
	{

		public:

			LinearElement2D GetElement(void);

		private:

			std::vector<std::vector<int>,std::vector<int>> GetStiffnessMatrix(void);

			void GetMassMatrix(void);
	};

}

