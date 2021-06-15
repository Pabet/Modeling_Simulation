//
// Created by Victor Stroescu on 06/11/2018.
//
#include "pvdWriter.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <list>
using namespace std;
namespace outputWriter {

    pvdWriter::pvdWriter(){}

    pvdWriter::~pvdWriter(){}

    void pvdWriter::writeFile(const std::string &file, std::list<std::string> record) {
        string outputanim=file;
        ofstream fid (outputanim);///@brief PVD file Generation based on all the VTU fles
        fid<<"<?xml version=\"1.0\"?>"<<endl;
        fid<<"<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">"<<endl;
        fid<<"<Collection>"<<endl;
        while (!record.empty())
        { std::string current =record.front();
        fid<<current;
        record.pop_front();}
        fid<<"</Collection>"<<endl;
        fid<<"</VTKFile>"<<endl;
        fid.close();
    }
}