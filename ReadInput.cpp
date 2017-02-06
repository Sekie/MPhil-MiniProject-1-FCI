#include <iostream>
#include <vector>
#include <string>
#include <fstream>

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

class InputObj
{
    public:
        void GetInputName();
        void Set();
    private:
        int aElectrons;
        int bElectrons;
        int aOrbitals;
        int bOrbitals;
        std::string InputName;
};

void InputObj::GetInputName()
{
    std::cout << "FCI: Enter input filename:" << std::endl;
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
}

int main()
{
    return 0;
}