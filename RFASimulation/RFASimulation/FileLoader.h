#pragma once
#include <string>
#include <iostream>
#include "TriangleMesh.h"

namespace RFASimulation
{
	static class FileLoader
	{

	public:

		TriangleMesh LoadMeshFromFilePath(string path);


	};

}