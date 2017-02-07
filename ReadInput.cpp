#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <map>
#include "ReadInput.h"

// class ReadInput(int &aElectrons, int &aOrbitals, int &bElectrons, int &bOrbitals)
// {
//     std::cout << "Enter input filename:" << std::endl;
//     std::string InputName;
//     std::cin >> InputName;
//     std::ifstream Input(InputName.c_str());

//     std::vector<double> tmpVec(4);
//     int 
//     while(Input >> tmpVec[0] >> tmpVec[1] >> tmpVec[2] >> tmpVec[3])
//     {
//         Molecule.push_back(tmpVec);
//     }
//     return 0;
// }

// class InputObj
// {
//     public:
//         void GetInputName();
//         void Set();
//     private:
//         int aElectrons;
//         int bElectrons;
//         int aOrbitals;
//         int bOrbitals;
//         std::string InputName;
//         std::map< std::string, double > Integrals;
// };

void InputObj::GetInputName()
{
    std::cout << "FCI: Enter input filename:\nFCI: ";
    std::cin >> InputName;
}

void InputObj::Set()
{
    std::ifstream InputFile(InputName.c_str());
    InputFile >> aElectrons >> aOrbitals >> bElectrons >> bOrbitals;
    if(aElectrons > aOrbitals || bElectrons > bOrbitals)
    {
        std::cerr << "FCI: ***** ERROR *****\nFCI: More electrons than orbitals. Check input file." << std::endl;
    }
    double tmpDouble;
    int tmpInt1, tmpInt2, tmpInt3, tmpInt4;
    while(!InputFile.eof())
    {
        /* We have to include all 8-fold permuation symmetries. This holds each integral in chemistry notation. We represent
        (ij|kl) as "i j k l". h_ij is "i j 0 0", as given in QChem. */
        InputFile >> tmpDouble >> tmpInt1 >> tmpInt2 >> tmpInt3 >> tmpInt4;
        Integrals[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
        Integrals[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
        Integrals[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
        Integrals[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
        Integrals[std::to_string(tmpInt2) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt4)] = tmpDouble;
        Integrals[std::to_string(tmpInt4) + " " + std::to_string(tmpInt3) + " " + std::to_string(tmpInt1) + " " + std::to_string(tmpInt2)] = tmpDouble;
        Integrals[std::to_string(tmpInt1) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt3)] = tmpDouble;
        Integrals[std::to_string(tmpInt3) + " " + std::to_string(tmpInt4) + " " + std::to_string(tmpInt2) + " " + std::to_string(tmpInt1)] = tmpDouble;
    }
}