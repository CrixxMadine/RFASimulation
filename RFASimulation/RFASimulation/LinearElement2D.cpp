#include "LinearElement2D.h"

using namespace RFASimulation; 



LinearElement2D::LinearElement2D(vector<double> k_coeff, vector<double> q_coeff, vector<double> f_rhs, 
	vector<vector<double>> coordinates, int numQuadPoints, bool buildMassMatrix)
{
	*m_p1 = coordinates[0].data();
	*m_p2 = coordinates[1].data();
	*m_p3 = coordinates[2].data();

	*m_kco = k_coeff.data();
	*m_qco = q_coeff.data();
	*m_frhs = f_rhs.data();

	generateJacobi();

	createStiffnessMatrix();
	createRightHandSide();

	if (buildMassMatrix) 
	{
		createMassMatrix();
	}
}

LinearElement2D::~LinearElement2D(void)
{
}

void LinearElement2D::generateJacobi()
{
	// Jacobi = [ (p2x - p1x), (p3x - p1x)
	//			  (p2y - p1y), (p3y - p1y) ]

	double a, b, c, d, det;

	a = *m_p2[0] - *m_p1[0];
	b = *m_p3[0] - *m_p1[0];
	c = *m_p2[1] - *m_p1[1];
	d = *m_p3[1] - *m_p1[1];

	det = a * d - c * b;


	// Inverted Jacobi =  1/det [ d  -b]
	//                          [-c   a]

	m_invJacobi[0][0] =  d / det;
	m_invJacobi[0][1] = -b / det;
	m_invJacobi[1][0] = -c / det;
	m_invJacobi[1][1] =  a / det;
}

void LinearElement2D::createStiffnessMatrix()
{
}

void LinearElement2D::createRightHandSide()
{
}


void LinearElement2D::createMassMatrix()
{
}
