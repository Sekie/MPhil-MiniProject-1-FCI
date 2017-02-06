#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <vector>
#include <time.h>
#include <utility> // Pair
#include <tuple>

void Davidson(Eigen::SparseMatrix<double> Ham, int Dim, int NumberOfEV, int L);

int BinomialCoeff(int n, int k) // n choose k
{
    int nCk = 1;
    int denom = 1;
    if (k <= 0 || k >= n)
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

/* Indexing equation given in Knowles and Handy */
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

/* This imposes an order onto the binary strings. The function takes an index and returns the corresponding binary string.
   We order the strings as such:
   0: 11000
   1: 10100
   2: 10010
   3: 10001
   4: 01100
   5: 01010
   6: 01001
   7: 00110
   8: 00101
   9: 00011
   We find the binary string recursively. We first ask, is the first digit 0 or 1. The first digit turns into a zero when
   all of the remaining electrons (n_electrons - 1) have permuted through all of the orbitals excluding the lowest one 
   (n_orbitals - 1). So if the index is above (n_orb - 1) choose (n_e - 1), then the first digit is zero. In this case, we
   ask ourselves the same problem for the remaining digits, but now we have one less orbital and we should subtract the number
   of permutations from the index. If the index is below, then the first digit is one. We fill that in and ask ourselves the 
   same question with the remaining digits. It is the same problem but with one less orbital AND one less electron, but the 
   index should not be changed. */
void GetOrbitalString(int Index, int NumElectrons, int NumOrbitals, std::vector<bool> &OrbitalString)
{
    if(NumOrbitals > 0) // Stop when we have chosen a digit for all orbitals. We take off an orbital each time we fill it.
    {
        int PreviousComb = BinomialCoeff(NumOrbitals - 1, NumElectrons - 1); // Number of ways for the higher electrons to permute in the higher orbitals.
        if (NumElectrons < 1) // If we don't have any electrons left, then all the remaining orbitals have to be empty.
        {
            OrbitalString.push_back(false);
            GetOrbitalString(Index, NumElectrons, NumOrbitals - 1, OrbitalString); // Move onto next orbital, remove one orbital from the list.
        }
        else if(Index < PreviousComb) // Means we have not finished all permutations and there is still an electron in the lowest orbital.
        {
            OrbitalString.push_back(true); // Put this electron there.
            GetOrbitalString(Index, NumElectrons - 1, NumOrbitals - 1, OrbitalString); // Consider the same problem, but with one less orbital and one less electron.
        }
        else // Means we have finished all those permuations and the electron in the first orbital should have moved.
        {
            Index -= PreviousComb; // Truncate the index, since we are considering a reduced problem.
            OrbitalString.push_back(false); // Empty orbital.
            GetOrbitalString(Index, NumElectrons, NumOrbitals - 1, OrbitalString); // Consider the same problem, but with one less orbital and a truncated index.
        }
    }
}

/* This counts the number of differences between two different binary representations */
int CountDifferences(std::vector<bool> BraString, std::vector<bool> KetString)
{
    int NumDiff = 0;
    for(int i = 0; i < KetString.size(); i++)
    {
        /* We want to count different orbitals, but if we end up doing that brute force we will double count since
           01 and 10 are different by one orbital, but their strings mismatch in two positions. It is enough to check
           that a 1 on one of the strings mismatches with the other string. */
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

/* 
   This is a function that finds the sign after putting two determinants into Slater-Condon form.
   To do so, we cound the number of positions each similar element is away from each other, and then
   we total this number. If the total is even, then the sign is positive and negative otherwise. Since
   the orbitals are organized ascending orbital number, we can always permute the lowest match, then the
   second lowest, and so on until they all line up. This justifies this method.
*/
int FindSign(std::vector<bool> BraString, std::vector<bool> KetString)
{
    int Sign = 0; // Will be +1 or -1, the sign.
    std::vector<int> KetOrder; // Counts position of the matching element in the ket.
    std::vector<int> BraOrder; // Does the same for the bra.
    int KetCount = 0;
    for(int i = 0; i < KetString.size(); i++)
    {
        if(KetString[i]) // We loop through both binary strings until we find a match, adding one to the counter each time to count our position.
        {
            KetCount++; // Count position along the ket. It only counts 1's, so filled orbitals.
            if(BraString[i]) // We have a match, but we don't know which position it is on the bra.
            {
                int BraCount = 1; // We go through the bra to count the number of orbitals (1's) before it. Not counting the actual element itself, so ammend with plus one.
                for(int j = 0; j < i; j++) // Loop through the string.
                {
                    if(BraString[j]) // Add one to the counter whereever an orbital is occupied.
                    {
                        BraCount++;
                    }
                }
                KetOrder.push_back(KetCount); // Add these positions to the list.
                BraOrder.push_back(BraCount);
            }
        }
    }

    for(int i = 0; i < KetOrder.size(); i++)
    {
        Sign += (abs(KetOrder[i] - BraOrder[i])); // Take hte total sum of position differences and use parity to find sign.
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

/* This converts a binary string into a string of orbital numbers. For example, 01011 -> 245 */
std::vector<int> ListOrbitals(std::vector<bool> DeterminantString)
{
    std::vector<int> OrbitalList;
    for(int i = 0; i < DeterminantString.size(); i++) // Loop through the binary string
    {
        if(DeterminantString[i]) // Put down the position of 1 for each 1 in the string.
        {
            OrbitalList.push_back(i + 1);
        }
    }
    return OrbitalList;
}

/* This function will take two binary strings and give the numerical number of the orbitals in which the two strings differ.
   For example, 0010110 and 0010011 differ in one orbital and can be reformulated as 356 and 367. We see that these two determinants
   disagree in orbital 5 and orbital 7. So when looking at <356|H|367>, this funciton stores 5, 7. That is, the difference from the 
   bra as the orbital number then the distant from the ket as the orbital number. For two electron differences, <1234|H|3456> is 
   listed in the vector as 1, 2, 5, 6, similar to the two electron integral denotation <ij|kl>. */ 
std::vector<int> ListDifference(std::vector<bool> BraString, std::vector<bool> KetString) // For 2e difference, put in <ij|kl> form
{
    std::vector<int> OrbitalList; // Holds the bra differences, we ammend the ket differences on afterwards.
    std::vector<int> tmpVec; // Holds the ket differences.
    for(int i = 0; i < KetString.size(); i++)
    {
        if(BraString[i] && !KetString[i]) // There is a mismatch, with the occupied orbital belonging to the bra. <..1..|..0..>
        {
            OrbitalList.push_back(i + 1);
        }
        if(!BraString[i] && KetString[i]) // Mismatch, occupied on the ket. <..0..|..1..>
        {
            tmpVec.push_back(i + 1);
        }
    }
    OrbitalList.insert(OrbitalList.end(), tmpVec.begin(), tmpVec.end()); // Stick the ket list onto the end.
    return OrbitalList;
}

int main()
{
    int aElectrons = 3;
    int bElectrons = 3;
    int aOrbitals = 10;
    int bOrbitals = 10;
    int NumberOfEV = 5;
    int L = NumberOfEV * 2;
    // int aVirtual = aOrbitals - aElectrons;
    // int bVirtual = bOrbitals - bElectrons;
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
        /* Since I order the orbitals the same way in the bra and ket, there should be no sign change */
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

    Start = clock();
    std::cout << "FCI: Beginning Direct Diagonalization... ";
    Eigen::MatrixXd HamDense = Ham;
    Eigen::SelfAdjointEigenSolver< Eigen::MatrixXd > HamEV;
    HamEV.compute(HamDense);
    std::cout << " done" << std::endl;
    std::cout << "FCI: The eigenvalues are\n" << HamEV.eigenvalues() << std::endl;
    std::cout << "FCI: Direct Diagonalization took " << (clock() - Start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    Start = clock();
    std::cout << "FCI: Beginning Davidson Diagonalization... " << std::endl;
    Davidson(Ham, Dim, NumberOfEV, L);
    std::cout << "...done" << std::endl;
    std::cout << "FCI: Davidson Diagonalization took " << (clock() - Start) / CLOCKS_PER_SEC << " seconds." << std::endl;

    return 0;
}