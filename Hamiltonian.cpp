#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <vector>

int BinomialCoeff(int n, int k)
{
    int nCk = 1;
    int denom = 1;
    if (k <= 0)
    {
        return 1;
    }
    for(int i = 0; i < k; i++)
    {
        nCk *= (n - i);
        denom *= (i + 1);
    }
    nCk /= denom;
    return nCk;
}

int Z_ForIndex(int ElectronNumber, int OrbitalNumber, int NumElectrons, int NumOrbitals)
{
    if(ElectronNumber == NumElectrons)
    {
        return OrbitalNumber - NumElectrons;
    }
    int Z = 0;
    for(int m = NumOrbitals - OrbitalNumber + 1; m <= NumOrbitals - ElectronNumber; m++)
    {
        Z += (BinomialCoeff(m, NumElectrons - ElectronNumber) - BinomialCoeff(m - 1, NumElectrons - ElectronNumber - 1));
    }
    return Z;
}

int StringIndex(std::vector<int> Orbitals, int NumOrbitals)
{
    int idx = 0;
    for (int i = 0; i < Orbitals.size(); i++)
    {
        idx += Z_ForIndex(i + 1, Orbitals[i], Orbitals.size(), NumOrbitals);
    }
    return idx;
}

void GetOrbitalString(int Index, int NumElectrons, int NumOrbitals, std::vector<bool> &OrbitalString)
{
    if(NumOrbitals > 0)
    {
        int PreviousComb = BinomialCoeff(NumOrbitals - 1, NumElectrons - 1);
        if (NumElectrons < 1)
        {
            OrbitalString.push_back(false);
            GetOrbitalString(Index, NumElectrons, NumOrbitals - 1, OrbitalString);
        }
        else if(Index < PreviousComb)
        {
            OrbitalString.push_back(true);
            GetOrbitalString(Index, NumElectrons - 1, NumOrbitals - 1, OrbitalString);
        }
        else
        {
            Index -= PreviousComb;
            OrbitalString.push_back(false);
            GetOrbitalString(Index, NumElectrons, NumOrbitals - 1, OrbitalString);
        }
    }
}

int CountDifferences(std::vector<bool> KetString, std::vector<bool> BraString)
{
    int NumDiff = 0;
    for(int i = 0; i < KetString.size(); i++)
    {
        if(KetString[i])
        {
            if(!BraString[i])
            {
                NumDiff++;
            }
        }
    }
    return NumDiff;
}

int FindSign(std::vector<bool> KetString, std::vector<bool> BraString)
{
    int Sign = 0;
    std::vector<int> KetOrder;
    std::vector<int> BraOrder;
    int KetCount = 0;
    for(int i = 0; i < KetString.size(); i++)
    {
        if(KetString[i])
        {
            KetCount++;
            if(BraString[i]) // Matching elements.
            {
                int BraCount = 1; // Not counting the actual element itself, so ammend with plus one.
                for(int j = 0; j < i; j++)
                {
                    if(BraString[j])
                    {
                        BraCount++;
                    }
                }
                KetOrder.push_back(KetCount);
                BraOrder.push_back(BraCount);
            }
        }
    }

    for(int i = 0; i < KetOrder.size(); i++)
    {
        Sign += (abs(KetOrder[i] - BraOrder[i]));
    }

    if(Sign % 2 == 0)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

int main()
{
    // OS.push_back(false);
    // foo(OS);
    // for(int i = 0; i < OS.size(); i++)
    // {
    //     std::cout << OS[i];
    // }
    // std::cout << std::endl;
    // return 0;

    int max = BinomialCoeff(7,4);

    std::vector<bool> OS1;
    std::vector<bool> OS2;

    for(int i = 0; i < max; i++) 
    {
        std::vector<bool> OS;
        int test = 0;
        GetOrbitalString(i, 4, 7, OS);
        std::cout << i << "\t";
        for(int j = 0; j < OS.size(); j++)
        {
            std::cout << OS[j];
        }
        std::cout << std::endl;

        if(i == 5)
        {
            OS1 = OS;
        }
        if(i == 20)
        {
            OS2 = OS;
        }
    }

    int Sign = FindSign(OS1, OS2);
    std::cout << Sign << std::endl;
    return 0;

    int aElectrons = 3;
    int bElectrons = 3;
    int aOrbitals = 15;
    int bOrbitals = 15;
    int aVirtual = aOrbitals - aElectrons;
    int bVirtual = bOrbitals - bElectrons;
    int aDim = BinomialCoeff(aOrbitals, aElectrons);
    int bDim = BinomialCoeff(bOrbitals, bElectrons);
    int Dim = aDim * bDim;

    std::vector< std::vector<bool> > aStrings;
    std::vector< std::vector<bool> > bStrings;
    for(int i = 0; i < aDim; i++)
    {
        std::vector<bool> tmpVec;
        GetOrbitalString(i, aElectrons, aOrbitals, tmpVec);
        aStrings.push_back(tmpVec);
    }
    for(int i = 0; i < bDim; i++)
    {
        std::vector<bool> tmpVec;
        GetOrbitalString(i, bElectrons, bOrbitals, tmpVec);
        bStrings.push_back(tmpVec);
    }

    Eigen::SparseMatrix<double> Ham(Dim, Dim);
    for(int i = 0; i < Dim; i++)
    {
        /* Diagonal Elements */
    }
    for(int i = 0; i < Dim; i++)
    {
        int aRow = i % bDim;
        int bRow = i / bDim;
        for(int j = i + 1; j < Dim; j++)
        {
            int aCol = j % bDim;
            int bCol = j / bDim; // Floor division
            if(CountDifferences(aStrings[aRow], aStrings[aCol]) + CountDifferences(bStrings[bRow], bStrings[bCol]) == 1)
            {

            }
            if(CountDifferences(aStrings[aRow], aStrings[aCol]) + CountDifferences(bStrings[bRow], bStrings[bCol]) == 2)
            {
                
            }
            else
            {
                
            }
        }
    }

    return 0;
}