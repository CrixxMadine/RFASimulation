#pragma once
#include <vector>


namespace RFASimulation
{
	using namespace std;


	class LinearElement2D
	{

		public: 
			
			vector<vector<int>> StiffnessMatrix;

			vector<vector<int>> MassMatrix;

			vector<int> RightHandSide;
			
			
			LinearElement2D(vector<double> k_coeff, vector<double> q_coeff, vector<double> f_rhs, 
				vector<vector<double>> coordinates, int numQuadPoints, bool buildMassMatrix);

			~LinearElement2D(void);

		private:

			double* m_kco[3];
			double* m_qco[3];
			double* m_frhs[3];

			double* m_p1[2];
			double* m_p2[2];
			double* m_p3[2];
			double m_invJacobi[2][2];

			static const int m_drPhi1 = -1;
			static const int m_drPhi2 =  1;
			static const int m_drPhi3 =  0;

			static const int m_dzPhi1 = -1;
			static const int m_dzPhi2 =  0;
			static const int m_dzPhi3 =  1;

			void generateJacobi();

			void createRightHandSide();

			void createStiffnessMatrix();

			void createMassMatrix();
	};

}

