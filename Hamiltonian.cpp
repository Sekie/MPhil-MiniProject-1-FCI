#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <vector>
#include <time.h>
#include <utility> // Pair

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

    clock_t Start = clock();

    std::cout << "FCI: Generating all determinant binary representations and enumerating determinants with differences... ";
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

    std::vector< std::pair<int, int> > aSingleDifference;
    std::vector< std::pair<int, int> > aDoubleDifference;
    int tmpInt;
    std::pair<int, int> tmpPair;
    for(int i = 0; i < aDim; i++)
    {
        for(int j = i + 1; j < aDim; j++)
        {
            tmpInt = CountDifferences(aStrings[i], aStrings[j]);
            if(tmpInt == 1)
            {
                tmpPair = std::make_pair(i,j);
                aSingleDifference.push_back(tmpPair);
            }
            if(tmpInt == 2)
            {
                tmpPair = std::make_pair(i,j);
                aDoubleDifference.push_back(tmpPair);

            }
        }
    }
    std::vector< std::pair<int, int> > bSingleDifference;
    std::vector< std::pair<int, int> > bDoubleDifference;
    for(int i = 0; i < bDim; i++)
    {
        for(int j = i + 1; j < bDim; j++)
        {
            tmpInt = CountDifferences(bStrings[i], bStrings[j]);
            if(tmpInt == 1)
            {
                tmpPair = std::make_pair(i,j);
                bSingleDifference.push_back(tmpPair);
            }
            if(tmpInt == 2)
            {
                tmpPair = std::make_pair(i,j);
                bDoubleDifference.push_back(tmpPair);

            }
        }
    }

    std::cout << "done.\nFCI: Commencing with matrix initialization... " << std::endl;;
    Eigen::SparseMatrix<double> Ham(Dim, Dim);

    typedef Eigen::Triplet<double> T;
    std::vector<T> tripletList;
    tripletList.reserve(Dim);

    // #pragma omp parallel for // Probably not worth diagonalizing since we have to include a critical in each loop.
    for(int i = 0; i < Dim; i++)
    {
        /* Diagonal Elements */
        tripletList.push_back(T(i, i , 1));
    }
    std::cout << "FCI: ...diagonal elements complete." << std::endl;

    /* 
       Now we begin setting the nonzero off-diagonal elements. We separate this into three groups.
          1) Elements differing by one spin-orbital (one orbital in alpha and beta electrons)
             This is achieved by looping through single differences of alpha (beta) and choosing
             to keep beta (alpha) diagonal so that there is only one difference.
          2) Elements differing by two spin-orbitals of the same spin.
             This is achieved by looping through double differences of alpha (beta) and choosing
             to keep beta (alpha) diagonal so that there is only two differences.
          3) Elements differing by two spin-orbitals of different spin.
             This is achieved by looping through single differences of both alpha and beta spaces
             and gives us two differences composed of a one different alpha and one different beta
             orbital.
    */

    /* 
       We start with Group 1. The matrix of single differences in alpha is block diagonal in beta states,
       but off diagonal in alpha states within that block.
               |      *   * |            |            |
               |        *   |            |            |
               |          * |            |            |
               |            |            |            |
               |____________|____________|____________|
               |            |      *   * |            |
               |            |        *   |            |
               |            |          * |            |
               |            |            |            |
               |____________|____________|____________|
               |            |            |     *   *  |
               |            |            |       *    |
               |            |            |         *  |
               |            |            |            |
               |____________|____________|____________|
    */
    #pragma omp parallel for
    for(int i = 0; i < aSingleDifference.size(); i++)
    {
        std::vector<T> tripletList_Private;
        int Index1, Index2;
        for(int j = 0; j < bDim; j++)
        {
            Index1 = aSingleDifference[i].first + j * aDim; // Diagonal in beta states. Hop to other beta blocks.
            Index2 = aSingleDifference[i].second + j * aDim;

            tripletList_Private.push_back(T(Index1, Index2 , 1));
            tripletList_Private.push_back(T(Index2, Index1 , 1));
            //std::cout << Index1 << "\t" << Index2 << std::endl;
        }
        #pragma omp critical
        tripletList.insert(tripletList.end(), tripletList_Private.begin(), tripletList_Private.end());
    }

    /* 
       The matrix of single differences in beta is not block diagonal, but is diagonal within the
       beta blocks that are off diagonal.
               |            |            |            |
               |            |  *         |            |
               |            |     *      |            |
               |            |        *   |            |
               |____________|____________|____________|
               |            |            |            |
               |            |            |  *         |
               |            |            |     *      |
               |            |            |        *   |
               |____________|____________|____________|
               |            |            |            |
               |            |            |            |
               |            |            |            |
               |            |            |            |
               |____________|____________|____________|
    */
    #pragma omp parallel for
    for(int i = 0; i < bSingleDifference.size(); i++)
    {
        std::vector<T> tripletList_Private;
        int Index1, Index2;
        for(int j = 0; j < aDim; j++)
        {
            Index1 = bSingleDifference[i].first * aDim + j; // Loop through each same alpha state in each beta block.
            Index2 = bSingleDifference[i].second * aDim + j;

            tripletList_Private.push_back(T(Index1, Index2 , 1));
            tripletList_Private.push_back(T(Index2, Index1 , 1));
        }
        #pragma omp critical
        tripletList.insert(tripletList.end(), tripletList_Private.begin(), tripletList_Private.end());
    }

    std::cout << "FCI: ...elements differing by one spin-orbital complete." << std::endl;

    /* Now Group 3. The elements of the matrix for two differences, exclusively alpha or beta spin-orbitals, has the same
       matrix form as before. We have to loop through the other spins having no differences. */
    #pragma omp parallel for
    for(int i = 0; i < aDoubleDifference.size(); i++)
    {
        std::vector<T> tripletList_Private;
        int Index1, Index2;
        for(int j = 0; j < bDim; j++)
        {
            Index1 = aDoubleDifference[i].first + j * aDim;
            Index2 = aDoubleDifference[i].second + j * aDim;

            tripletList_Private.push_back(T(Index1, Index2 , 1));
            tripletList_Private.push_back(T(Index2, Index1 , 1));
        }
        #pragma omp critical
        tripletList.insert(tripletList.end(), tripletList_Private.begin(), tripletList_Private.end());
    }
    #pragma omp parallel for
    for(int i = 0; i < bDoubleDifference.size(); i++)
    {
        std::vector<T> tripletList_Private;
        int Index1, Index2;

        for(int j = 0; j < aDim; j++)
        {
            Index1 = bDoubleDifference[i].first * aDim + j; // Loop through each same alpha state in each beta block.
            Index2 = bDoubleDifference[i].second * aDim + j;

            tripletList_Private.push_back(T(Index1, Index2 , 1));
            tripletList_Private.push_back(T(Index2, Index1 , 1));
        }
        #pragma omp critical
        tripletList.insert(tripletList.end(), tripletList_Private.begin(), tripletList_Private.end());
    }

    /* Now Group 3. Unlike before, we don't have to loop over alpha or beta having no differences. We simply loop
       over both alpha and beta having one difference. */
    #pragma omp parallel for
    for(int i = 0; i < aSingleDifference.size(); i++)
    {
        std::vector<T> tripletList_Private;
        int Index1, Index2;
        for(int j = 0; j < bSingleDifference.size(); j++)
        {
            Index1 = aSingleDifference[i].first + aDim * bSingleDifference[j].first;
            Index2 = aSingleDifference[i].second + aDim * bSingleDifference[j].second;
            //Ham.insert(Index1, Index2) = 1;
            //Ham.insert(Index2, Index1) = 1;
            tripletList_Private.push_back(T(Index1, Index2 , 1));
            tripletList_Private.push_back(T(Index2, Index1 , 1));
            /* We have to be a little more careful in this case. We want the upper triangle, but this only gives us half 
               of the upper triangle. In particular, the upper half of each beta block in upper triangle of the full matrix
               are the only nonzero elements. We want the whole beta block in the upper triangle of the full matrix to be
               nonzero where needed.
               |            |      *   * |      *   * |
               |            |        *   |        *   |
               |            |          * |          * |
               |            |            |            |
               |____________|____________|____________|
               |            |            |      *   * |
               |            |            |        *   |
               |            |            |          * |
               |            |            |            |
               |____________|____________|____________|
               |            |            |            |
               |            |            |            |
               |            |            |            |
               |            |            |            |
               |____________|____________|____________|
               So we have to include some transposed elements too. It is enough to transpose the alpha indices. this
               transposes each block above, and we end up with a fully upper triangular matrix. */
            Index1 = aSingleDifference[i].second + aDim * bSingleDifference[j].first;
            Index2 = aSingleDifference[i].first + aDim * bSingleDifference[j].second; // Note that first and second are switched for alpha here.

            tripletList_Private.push_back(T(Index1, Index2 , 1));
            tripletList_Private.push_back(T(Index2, Index1 , 1));
            // std::cout << Index1 << "\t" << Index2 << std::endl;
        }
        #pragma omp critical
        tripletList.insert(tripletList.end(), tripletList_Private.begin(), tripletList_Private.end());
    }
    
    std::cout << "FCI: ...elements differing by two spin-orbitals complete." << std::endl;

    Ham.setFromTriplets(tripletList.begin(), tripletList.end());

    std::cout << "FCI: Hamiltonian initialization took " << (clock() - Start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return 0;
}