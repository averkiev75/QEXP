// QEXP.cpp : Version 2 1/27/2017
// Now takes into account Card's options:	alat | bohr | angstrom | crystal 

#include "stdafx.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <regex>
#include <map>
#define _USE_MATH_DEFINES
#include <math.h>
using namespace std;

string line;
fstream infile, outfile;
float x, y, z;
float a, b, c, al, bt, gm, det;
float matr[3][3], invt[3][3];
float(*coord)[3];
string  *aname;
string  *element;
vector<string> types;
const float Bohr = 0.529177249;
float cell(float(*matr)[3], float(*invt)[3], float *a, float *b, float *c, float *al, float *bt, float *gm);
void makeOutput(string filename, int structCurrent);
int atoms = 0; //total number of atoms
int scale = 1;
float alat = 1;

int main(int argc, char **argv)
{
	string filename;
	int structCurrent = 0, structNumber = 9999;
	char *next;
	if (argc == 1) {
		cout << "Usage: QEXP filename - extract the last structure or\n";
		cout << "QEXP filename N - extract structure from step N\n";
		cout << "for N = 0, program extracts initial structure\n";
		cout << "for negative N, program extracts all structures\n";
		return 1;
	}
	filename = argv[1];
	if (argc > 2)
		structNumber = strtof(argv[2], NULL);
	types.push_back("C");
	types.push_back("H");
	infile.open(filename);
	if (infile.is_open()) {
		do { getline(infile, line); } while (int(line.find(" lattice parameter (alat)")) < 0);
		alat = strtof(line.c_str() + 34, NULL); alat *= Bohr;
		do { getline(infile, line); } while (int(line.find("number of atoms")) < 0);
		atoms = strtof(line.c_str() + 34, NULL);
		coord = new float[atoms][3];
		aname = new string[atoms];
		element = new string[atoms];
		do { getline(infile, line); } while (int(line.find("crystal axes: (cart. coord. in units of alat)")) < 0);
		getline(infile, line);
		matr[0][0] = alat*strtof(line.c_str() + 24, &next);
		matr[0][1] = alat*strtof(next, &next);
		matr[0][2] = alat*strtof(next, NULL);
		getline(infile, line);
		matr[1][0] = alat*strtof(line.c_str() + 24, &next);
		matr[1][1] = alat*strtof(next, &next);
		matr[1][2] = alat*strtof(next, NULL);
		getline(infile, line);
		matr[2][0] = alat*strtof(line.c_str() + 24, &next);
		matr[2][1] = alat*strtof(next, &next);
		matr[2][2] = alat*strtof(next, NULL);
		det = cell(matr, invt, &a, &b, &c, &al, &bt, &gm);
		printf(" 0 %6.4f  %6.4f  %6.4f  %6.3f  %6.3f  %6.3f  %6.3f\n", a, b, c, al, bt, gm, det);
		do { getline(infile, line); } while (int(line.find(" site n.     atom                  positions (alat units)")) < 0);
		getline(infile, line);
		for (int at = 0; at < atoms; at++) {
			aname[at] = string(line, 21, 4);
			coord[at][0] = strtof(line.c_str() + 40, &next);
			coord[at][1] = strtof(next, &next);
			coord[at][2] = strtof(next, NULL);
			getline(infile, line);
		}
		if (structNumber != 0) {
			if (structNumber < 0)
				makeOutput(filename, 0);
			while (getline(infile, line)) {
				if (int(line.find("CELL_PARAMETERS")) > -1) {
					getline(infile, line);
					matr[0][0] = strtof(line.c_str(), &next);
					matr[0][1] = strtof(next, &next);
					matr[0][2] = strtof(next, NULL);
					getline(infile, line);
					matr[1][0] = strtof(line.c_str(), &next);
					matr[1][1] = strtof(next, &next);
					matr[1][2] = strtof(next, NULL);
					getline(infile, line);
					matr[2][0] = strtof(line.c_str(), &next);
					matr[2][1] = strtof(next, &next);
					matr[2][2] = strtof(next, NULL);
					det = cell(matr, invt, &a, &b, &c, &al, &bt, &gm);
					printf("%2i %6.4f  %6.4f  %6.4f  %6.3f  %6.3f  %6.3f  %6.3f\n", structCurrent + 1, a, b, c, al, bt, gm, det);
				}
				else if (int(line.find("ATOMIC_POSITIONS")) > -1) {
					if (int(line.find("alat")) > -1)
						scale = 1;
					else if (int(line.find("bohr")) > -1)
						scale = 2;
					else if (int(line.find("angstrom")) > -1)
						scale = 3;
					else if (int(line.find("crystal")) > -1)
						scale = 4;
					else
						cout << "\nCheck ATOMIC_POSITIONS line\n";
					getline(infile, line);
					for (int at = 0; at < atoms; at++) {
						aname[at] = string(line, 0, 8);
						coord[at][0] = strtof(line.c_str() + 4, &next);
						coord[at][1] = strtof(next, &next);
						coord[at][2] = strtof(next, NULL);
						getline(infile, line);
					}
					structCurrent++;
					if (structCurrent == structNumber)
						break;
					else if (structNumber < 0)
						makeOutput(filename, structCurrent);
				}
			}
			makeOutput(filename, structCurrent);
		}
		else
			makeOutput(filename, 0);
	}
	else
		cout << "Cannot open file: " << argv[1] << endl;
	delete[] coord;
	delete[] aname;
	infile.close();
	cout << "Read " << structCurrent << " structures\n";
	return 0;
}


float cell(float(*matr)[3], float(*invt)[3], float *a, float *b, float *c, float *al, float *bt, float *gm)
{
	*a = sqrt(matr[0][0] * matr[0][0] + matr[0][1] * matr[0][1] + matr[0][2] * matr[0][2]);
	*b = sqrt(matr[1][0] * matr[1][0] + matr[1][1] * matr[1][1] + matr[1][2] * matr[1][2]);
	*c = sqrt(matr[2][0] * matr[2][0] + matr[2][1] * matr[2][1] + matr[2][2] * matr[2][2]);
	*al = acos((matr[1][0] * matr[2][0] + matr[1][1] * matr[2][1] + matr[1][2] * matr[2][2]) / (*b) / (*c)) * 180 / M_PI;
	*bt = acos((matr[0][0] * matr[2][0] + matr[0][1] * matr[2][1] + matr[0][2] * matr[2][2]) / (*a) / (*c)) * 180 / M_PI;
	*gm = acos((matr[1][0] * matr[0][0] + matr[1][1] * matr[0][1] + matr[1][2] * matr[0][2]) / (*b) / (*a)) * 180 / M_PI;
	float det = matr[0][0] * matr[1][1] * matr[2][2] + matr[0][1] * matr[1][2] * matr[2][0] + matr[0][2] * matr[1][0] * matr[2][1] -
		matr[0][0] * matr[1][2] * matr[2][1] - matr[0][1] * matr[1][0] * matr[2][2] - matr[0][2] * matr[1][1] * matr[2][0];
	invt[0][0] = (matr[1][1] * matr[2][2] - matr[1][2] * matr[2][1]) / det;
	invt[0][1] = -(matr[1][0] * matr[2][2] - matr[1][2] * matr[2][0]) / det;
	invt[0][2] = (matr[1][0] * matr[2][1] - matr[1][1] * matr[2][0]) / det;
	invt[1][0] = -(matr[0][1] * matr[2][2] - matr[0][2] * matr[2][1]) / det;
	invt[1][1] = (matr[0][0] * matr[2][2] - matr[0][2] * matr[2][0]) / det;
	invt[1][2] = -(matr[0][0] * matr[2][1] - matr[0][1] * matr[2][0]) / det;
	invt[2][0] = (matr[0][1] * matr[1][2] - matr[0][2] * matr[1][1]) / det;
	invt[2][1] = -(matr[0][0] * matr[1][2] - matr[0][2] * matr[1][0]) / det;
	invt[2][2] = (matr[0][0] * matr[1][1] - matr[0][1] * matr[1][0]) / det;
	return det;
}

void makeOutput(string filename, int structCurrent)
{
	for (int at = 0; at < atoms; at++) {
		element[at] = string(aname[at], 0, 2);
		element[at][0] = toupper(element[at][0]);
		if (isalpha(element[at][1])) {
			element[at][1] = tolower(element[at][1]);
		}
		else
			element[at] = string(aname[at], 0, 1);
		bool found = false;
		for (int i = 0; i < types.size(); i++) {
			if (types[i] == element[at]) {
				found = true;
				break;
			}
		}
		if (!found)
			types.push_back(element[at]);
	}

	int size = filename.size();
	filename = string(filename, 0, size - 4);
	if (structCurrent < 10) filename += "0";
	filename += to_string(static_cast<long long>(structCurrent));
	filename += ".res";
	outfile.open(filename, fstream::out);
	outfile << "TITL \nCELL 0.71073  " << fixed << setprecision(4) << a << "  " << b << "  " << c << "  " << setprecision(3) << al << "  " << bt << "  " << gm << endl;
	outfile << "LATT -1\nSFAC ";
	for (int i = 0; i < types.size(); i++)
		outfile << types[i] << " ";

	outfile << endl << endl;
	for (int at = 0; at < atoms; at++) {
		int type = 1;
		for (int i = 0; i < types.size(); i++) {
			if (types[i] == element[at]) {
				type = i + 1;
				break;
			}
		}
		if (scale != 4) {
			x = invt[0][0] * coord[at][0] + invt[0][1] * coord[at][1] + invt[0][2] * coord[at][2];
			y = invt[1][0] * coord[at][0] + invt[1][1] * coord[at][1] + invt[1][2] * coord[at][2];
			z = invt[2][0] * coord[at][0] + invt[2][1] * coord[at][1] + invt[2][2] * coord[at][2];
			if (scale == 1) { x *= alat; y *= alat; z *= alat; }
			else if (scale == 2) { x *= Bohr; y *= Bohr; z *= Bohr; }
		}
		else {
			x = coord[at][0];
			y = coord[at][1];
			z = coord[at][2];
		}

		element[at] += to_string(static_cast<long long>(at + 1));

		outfile << setw(7) << left << element[at] << type;
		outfile << right << fixed << setprecision(5) << setw(11) << x << setw(11) << y << setw(11) << z << "   11.00000    0.05\n";
	}
	outfile << "\nHKLF 4\nEND\n";
	outfile.close();
}
