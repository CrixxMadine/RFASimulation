#include "FileLoader.h"
#include <string>
#include <iostream>
#include "fstream"  
#include "sstream" 

using namespace std;

RFASimulation::TriangleMesh RFASimulation::FileLoader::LoadMeshFromFilePath(string path)
{
    return TriangleMesh();

    //int x, y;

    //ifstream inStream(path, ifstream::in);

    //if (!inStream) {
    //    cout << "Cannot open file.\n";
    //    return;
    //}

    //for (y = 0; y < 15; y++) {
    //    for (x = 0; x < 15; x++) {
    //        inStream >> distances[x][y];
    //    }
    //}

    //inStream.close();
}

static vector<vector<double>> LoadMatrixFromFile(string &filename)
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
            reader << val;

            if (reader.fail())
                break;

            lineData.push_back(val);
        }

        data.push_back(lineData);
    }
}

