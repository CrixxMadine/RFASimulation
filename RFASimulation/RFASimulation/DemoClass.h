#pragma once


namespace RFASimulation
{

	class DemoClass
	{
		public:

			DemoClass(int myInt);

			int GetMyInt(void);

			~DemoClass(void);

		private:

			int m_myInt = 0;

	};
}