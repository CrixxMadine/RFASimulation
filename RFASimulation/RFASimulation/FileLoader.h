#pragma once
#include <string>
#include <iostream>
#include "TriangleMesh.h"

namespace RFASimulation
{
	using namespace std;

	static class FileLoader
	{

	public:

		static TriangleMesh LoadMeshFromFilePath(string path);


	private:

		static vector<vector<double>> LoadMatrixFromFile(string filename);
		
	};

}