#include "FileLoader.h"
#include <string>
#include <iostream>
#include "fstream"  
#include "sstream" 

using namespace std;


RFASimulation::TriangleMesh RFASimulation::FileLoader::LoadMeshFromFilePath(string folderPath)
{
    string fileName = ("pmesh.txt");
    string fullPath = folderPath + fileName;

    auto pmesh = FileLoader::LoadMatrixFromFile(fullPath);

    fileName = ("tmesh.txt");
    fullPath = folderPath + fileName;

    auto tmesh_double = FileLoader::LoadMatrixFromFile(fullPath);

    // Convert double 2dvector to integer 2dvector
    vector<vector<int>> tmesh_int;
    tmesh_int.reserve(tmesh_double.size());

    for (auto&& v : tmesh_double)
    {
        tmesh_int.emplace_back(std::begin(v), std::end(v));
    }

    fileName = ("bedges.txt");
    fullPath = folderPath + fileName;

    auto bedges = FileLoader::LoadMatrixFromFile(fullPath);

    return TriangleMesh(tmesh_int, pmesh, bedges);

}

vector<vector<double>> RFASimulation::FileLoader::LoadMatrixFromFile(string filename)
{
    ifstream file;
    file.open(filename, ios::in | ios::out);

    if (!file.is_open()) {
        // error
        return vector<vector<double> >();
    }

    vector<vector<double> > data;
    string line;

    while (!std::getline(file, line, '\n').eof())
    {
        stringstream reader(line);

        vector<double> lineData;



        string::const_iterator i = line.begin();

        while (!reader.eof()) {

            double val;

            reader >> val;

            if (reader.fail())
                break;

            lineData.push_back(val);
        }

        data.push_back(lineData);
    }

    return data;
}

