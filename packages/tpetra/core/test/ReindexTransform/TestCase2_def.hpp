// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_2_DEF_HPP
#define TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_2_DEF_HPP

#include "TestCase2_decl.hpp"
#include <MatrixMarket_Tpetra.hpp>

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t>::TestCase2( Teuchos::RCP<Teuchos::Comm<int> const> comm )
  : TestCaseBase<Scalar_t, LocalId_t, GlobalId_t, Node_t>( comm
							 , 17 // globalNumRows
                                                         , 17 // globalNumCols
                                                         , 0  // rowIndexBase
                                                         , 0  // colIndexBase
                                                         , 5  // maxNnzPerRow
                                                         , 43 // globalNumNnz
                                                         , 4.940647893865e+00                       // frobNorm
                                                         , {0., 0.} // {1.000000073665e+00, 1.000000073665e+00} // lhsNorms2
                                                         , {0., 0.} // {1.000000075938e+00, 1.000000075938e+00} // rhsNorms2
                                                         )
{
  std::cout << "Entering TestCase2<>::constructor(), comm->getSize() = " << comm->getSize() << std::endl;
  if (comm->getSize() == 2) {
    std::vector< std::vector<GlobalId_t> > globalRowIndicesPerRank_beforeTransform;
    std::vector< std::vector<GlobalId_t> > globalColIndicesPerRank_beforeTransform;
    std::vector< std::vector<GlobalId_t> > globalDomainIndicesPerRank_beforeTransform;
    std::vector< std::vector<GlobalId_t> > globalRangeIndicesPerRank_beforeTransform;
    std::vector< std::vector< std::vector< ColValPair_t > > > colValuePairsPerRankAndLocalRow;
    this->prepareDataToCreateMatrix( globalRowIndicesPerRank_beforeTransform
                                   , globalColIndicesPerRank_beforeTransform
                                   , globalDomainIndicesPerRank_beforeTransform
                                   , globalRangeIndicesPerRank_beforeTransform
                                   , colValuePairsPerRankAndLocalRow
                                   );

    size_t globalNumNnz(0);
    for (int iRank(0); iRank < this->m_numRanks; ++iRank) {
      for (size_t j(0); j < colValuePairsPerRankAndLocalRow[iRank].size(); ++j) {
        globalNumNnz += colValuePairsPerRankAndLocalRow[iRank][j].size();
      }
    }
    if (globalNumNnz != this->m_globalNumNnz) {
      std::stringstream msg;
      msg << "In TestCase2::constructor()"
          << ": globalNumNnz should be " << this->m_globalNumNnz
          << ", but it is " << globalNumNnz
          << std::endl;
      throw std::runtime_error( msg.str() );
    }
    //std::cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++01" << std::endl;

    this->m_rowMap = std::unique_ptr<Map_t>( new Map_t( this->m_globalNumRows
                                                      , globalRowIndicesPerRank_beforeTransform[this->m_myRank].data()
                                                      , globalRowIndicesPerRank_beforeTransform[this->m_myRank].size()
                                                      , this->m_rowIndexBase
                                                      , this->m_comm
                                                      ));

    this->m_colMap = std::unique_ptr<Map_t>( new Map_t( this->m_globalNumCols
                                                      , globalColIndicesPerRank_beforeTransform[this->m_myRank].data()
                                                      , globalColIndicesPerRank_beforeTransform[this->m_myRank].size()
                                                      , this->m_colIndexBase
                                                      , this->m_comm
                                                      ));

    this->m_matrix = std::unique_ptr<Matrix_t>( new Matrix_t( Teuchos::rcp<Map_t>(this->m_rowMap.get(), false)
                                                            , Teuchos::rcp<Map_t>(this->m_colMap.get(), false)
                                                            , this->m_maxNnzPerRow
                                                            ));

    std::vector<Scalar_t>   tpetraValues (this->m_maxNnzPerRow);
    std::vector<GlobalId_t> tpetraIndices(this->m_maxNnzPerRow);
#if 0
    this->m_comm->barrier();

    std::cout.flush();
    this->m_comm->barrier();
    for (int p(0); p < this->m_numRanks; ++p) {
      this->m_comm->barrier();
      if (p != this->m_myRank) continue;

      std::cout.flush();
#endif

    for (size_t i(0); i < globalRowIndicesPerRank_beforeTransform[this->m_myRank].size(); ++i) {
      int globalRow = globalRowIndicesPerRank_beforeTransform[this->m_myRank][i];
      size_t numEntries = colValuePairsPerRankAndLocalRow[this->m_myRank][i].size();
      //std::cout << "p=" << p << ", calling m_matrix->insertGlobalValues() for i = " << i << ", numEntries = " << numEntries << std::endl;
      for (size_t j(0); j < numEntries; ++j) {
        tpetraValues [j] = colValuePairsPerRankAndLocalRow[this->m_myRank][i][j].second;
        tpetraIndices[j] = globalColIndicesPerRank_beforeTransform[this->m_myRank][ colValuePairsPerRankAndLocalRow[this->m_myRank][i][j].first ];
        //std::cout << "p,i,j,val,index=" << p << ", " << i << ", " << j << ", " << tpetraValues[j] << ", " << tpetraIndices[j] << std::endl;
      }

      this->m_matrix->insertGlobalValues( globalRow
                                        , numEntries
                                        , tpetraValues.data()
                                        , tpetraIndices.data()
                                        );
      //std::cout << "p=" << p << ", returned from m_matrix->insertGlobalValues() for i = " << i << std::endl;
    }

#if 0
    std::cout.flush();

    } // for p
    this->m_comm->barrier();
#endif

    this->m_matrix->fillComplete();
    this->m_matrixRCP = Teuchos::rcp<Matrix_t>(this->m_matrix.get(), false);

    size_t epetraNumVectors = 2;

    this->m_lhs = std::unique_ptr<MultiVector_t>( new MultiVector_t( Teuchos::rcp<Map_t>(this->m_rowMap.get(), false)
                                                                   , epetraNumVectors
                                                                   , true /* zeroOut */
                                                                   ));
    this->m_lhsRCP = Teuchos::rcp<MultiVector_t>(this->m_lhs.get(), false);

    this->m_rhs = std::unique_ptr<MultiVector_t>( new MultiVector_t( Teuchos::rcp<Map_t>(this->m_rowMap.get(), false)
                                                                   , epetraNumVectors
                                                                   , true /* zeroOut */
                                                                   ));
    this->m_rhsRCP = Teuchos::rcp<MultiVector_t>(this->m_rhs.get(), false);

#if 0
    Teuchos::ArrayRCP<scalar_t> tmpView = m_tpetraMultiVector->get1dViewNonConst();
    for (size_t v(0); v < epetraNumVectors; ++v) {
      double * const & myValues = epetraMultiVector[v];
      for (size_t i(0); i < epetraMultiVector.MyLength(); ++i) {
        tmpView[v * epetraMultiVector.MyLength() + i] = myValues[i];
      }
    }
#endif

    // Create linear problem
    this->m_linearProblem = std::unique_ptr<Problem_t>(new Problem_t(this->m_matrixRCP, this->m_lhsRCP, this->m_rhsRCP));

    this->m_matrixMapsShallChange = true;
  }
  else {
    this->baseInstantiateLinearProblem( "matA1.txt"
                                      , "lhs1.txt"
                                      , "rhs1.txt"
                                      );
    this->m_matrixMapsShallChange = false;
  }
  std::cout << "In TestCase2<>::constructor()" << std::endl;

  this->baseCheckBeforeOrAfterTransform( this->m_linearProblem.get() // Input
                                       , "Case2_original"            // Input
                                       );

  std::cout << "Leaving TestCase2<>::constructor()" << std::endl;
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t>::~TestCase2()
{
  // Nothing to do
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
bool
TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t>::checkTransformedProblem( Problem_t const * transformedProblem ) const
{
  return this->baseCheckTransformedProblem( transformedProblem );
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
bool
TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t>::checkAfterFwdRvs( Problem_t const * originalProblem ) const
{
  return this->baseCheckAfterFwdRvs( originalProblem );
}

template< class Scalar_t, class LocalId_t, class GlobalId_t, class Node_t > 
void
TestCase2<Scalar_t, LocalId_t, GlobalId_t, Node_t>::prepareDataToCreateMatrix
  ( std::vector< std::vector<GlobalId_t> >                  & globalRowIndicesPerRank
  , std::vector< std::vector<GlobalId_t> >                  & globalColIndicesPerRank
  , std::vector< std::vector<GlobalId_t> >                  & globalDomainIndicesPerRank
  , std::vector< std::vector<GlobalId_t> >                  & globalRangeIndicesPerRank
  , std::vector< std::vector< std::vector<ColValPair_t> > > & colValuePairsPerRankAndLocalRow
  )
{
  globalRowIndicesPerRank.resize(this->m_numRanks);
  globalColIndicesPerRank.resize(this->m_numRanks);
  colValuePairsPerRankAndLocalRow.resize(this->m_numRanks);

  if (this->m_numRanks == 2) {
    //RowMap    0 1 2 3 4 5 6 7 9
    //ColMap    0 1 2 3 4 5 6 7 9 12 13 11
    //DomainMap 0 1 2 3 4 5 6 7 9
    //RangeMap  0 1 2 3 4 5 6 7 9

    //RowMap    11 12 13 14 15 16 17 18
    //ColMap    11 12 13 14 15 16 17 18 9 7
    //DomainMap 11 12 13 14 15 16 17 18
    //RangeMap  11 12 13 14 15 16 17 18

    globalRowIndicesPerRank[0] = {0, 1, 2, 3, 4, 5, 6, 7, 9};
    globalRowIndicesPerRank[1] = {11, 12, 13, 14, 15, 16, 17, 18};

    globalColIndicesPerRank[0] = {0, 1, 2, 3, 4, 5, 6, 7, 9, 12, 13, 11};
    globalColIndicesPerRank[1] = {11, 12, 13, 14, 15, 16, 17, 18, 9, 7};

    for (int iRank(0); iRank < this->m_numRanks; ++iRank) {
      colValuePairsPerRankAndLocalRow[iRank].resize( globalRowIndicesPerRank[iRank].size() );
    }

    //Rank 0
    //i = 0, numEntries = 2; indices = 1 2; (0,1)=1. (0,2)=1.
    //i = 1, numEntries = 2; indices = 0 3; (1,0)=1. (1,3)=-1.
    //i = 2, numEntries = 2; indices = 0 4; (2,0)=1. (2,4)=-1.
    //i = 3, numEntries = 2; indices = 1 5; (3,1)=-1. (3,5)=-1.
    //i = 4, numEntries = 2; indices = 2 6; (4,2)=-1. (4,6)=1.
    //i = 5, numEntries = 2; indices = 3 7; (5,3)=-1. (5,7)=1.
    //i = 6, numEntries = 1; indices = 4; (6,4)=1.
    //i = 7, numEntries = 3; indices = 5 7 9; (7,5)=1. (7,7)=0.25 (7,9)=-0.25
    //i = 8, numEntries = 1; indices = 11; (8,11)=1.

    colValuePairsPerRankAndLocalRow[0][0].resize(2); // Rank 0 has 9 rows: local row 0 has 2 entries
    colValuePairsPerRankAndLocalRow[0][0][0] = ColValPair_t({1,1.});
    colValuePairsPerRankAndLocalRow[0][0][1] = ColValPair_t({2,1.});

    colValuePairsPerRankAndLocalRow[0][1].resize(2); // Rank 0 has 9 rows: local row 1 has 2 entries
    colValuePairsPerRankAndLocalRow[0][1][0] = ColValPair_t({0,1.});
    colValuePairsPerRankAndLocalRow[0][1][1] = ColValPair_t({3,-1.});

    colValuePairsPerRankAndLocalRow[0][2].resize(2); // Rank 0 has 9 rows: local row 2 has 2 entries
    colValuePairsPerRankAndLocalRow[0][2][0] = ColValPair_t({0,1.});
    colValuePairsPerRankAndLocalRow[0][2][1] = ColValPair_t({4,-1.});

    colValuePairsPerRankAndLocalRow[0][3].resize(2); // Rank 0 has 9 rows: local row 3 has 2 entries
    colValuePairsPerRankAndLocalRow[0][3][0] = ColValPair_t({1,-1.});
    colValuePairsPerRankAndLocalRow[0][3][1] = ColValPair_t({5,-1.});

    colValuePairsPerRankAndLocalRow[0][4].resize(2); // Rank 0 has 9 rows: local row 4 has 2 entries
    colValuePairsPerRankAndLocalRow[0][4][0] = ColValPair_t({2,-1.});
    colValuePairsPerRankAndLocalRow[0][4][1] = ColValPair_t({6,1.});

    colValuePairsPerRankAndLocalRow[0][5].resize(2); // Rank 0 has 9 rows: local row 5 has 2 entries
    colValuePairsPerRankAndLocalRow[0][5][0] = ColValPair_t({3,-1.});
    colValuePairsPerRankAndLocalRow[0][5][1] = ColValPair_t({7,1.});

    colValuePairsPerRankAndLocalRow[0][6].resize(1); // Rank 0 has 9 rows: local row 6 has 1 entry
    colValuePairsPerRankAndLocalRow[0][6][0] = ColValPair_t({4,1.});

    colValuePairsPerRankAndLocalRow[0][7].resize(3); // Rank 0 has 9 rows: local row 7 has 3 entries
    colValuePairsPerRankAndLocalRow[0][7][0] = ColValPair_t({5,1.});
    colValuePairsPerRankAndLocalRow[0][7][1] = ColValPair_t({7,0.25});
    colValuePairsPerRankAndLocalRow[0][7][2] = ColValPair_t({9,-0.25});

    colValuePairsPerRankAndLocalRow[0][8].resize(1); // Rank 0 has 9 rows: local row 8 has 1 entry
    colValuePairsPerRankAndLocalRow[0][8][0] = ColValPair_t({11,1.});

    //Rank 1
    //i = 0, numEntries = 4; indices = 0 1 5 8; (0,0)=9.29367e-05 (0,1)=-8.77169e-05 (0,5)=-5.21976e-06 (0,8)=1
    //i = 1, numEntries = 5; indices = 0 1 2 5 9; (1,0)=-5.6556e-35 (1,1)=0.25 (1,2)=-2.89557e-35 (1,5)=-3.29088e-35 (1,9)=-0.25
    //i = 2, numEntries = 3; indices = 1 2 4; (2,1)=0.0012461 (2,2)=0.000160258 (2,4)=-1
    //i = 3, numEntries = 2; indices = 3 6; (3,3)=0.2 (3,6)=-0.2
    //i = 4, numEntries = 2; indices = 2 6; (4,2)=-1 (4,6)=1
    //i = 5, numEntries = 4; indices = 0 1 5 7; (5,0)=-9.29367e-05 (5,1)=8.77169e-05 (5,5)=5.21976e-06 (5,7)=-1
    //i = 6, numEntries = 4; indices = 3 4 6 7; (6,3)=-0.2 (6,4)=1 (6,6)=0.2 (6,7)=1
    //i = 7, numEntries = 2; indices = 5 6; (7,5)=-1 (7,6)=1

    colValuePairsPerRankAndLocalRow[1][0].resize(4); // Rank 1 has 8 rows: local row 0 has 4 entries
    colValuePairsPerRankAndLocalRow[1][0][0] = ColValPair_t({0,9.29367e-05});
    colValuePairsPerRankAndLocalRow[1][0][1] = ColValPair_t({1,-8.77169e-05});
    colValuePairsPerRankAndLocalRow[1][0][2] = ColValPair_t({5,-5.21976e-06});
    colValuePairsPerRankAndLocalRow[1][0][3] = ColValPair_t({8,1.});

    colValuePairsPerRankAndLocalRow[1][1].resize(5); // Rank 1 has 8 rows: local row 1 has 5 entries
    colValuePairsPerRankAndLocalRow[1][1][0] = ColValPair_t({0,-5.6556e-35});
    colValuePairsPerRankAndLocalRow[1][1][1] = ColValPair_t({1,0.25});
    colValuePairsPerRankAndLocalRow[1][1][2] = ColValPair_t({2,-2.89557e-35});
    colValuePairsPerRankAndLocalRow[1][1][3] = ColValPair_t({5,-3.29088e-35});
    colValuePairsPerRankAndLocalRow[1][1][4] = ColValPair_t({9,-0.25});

    colValuePairsPerRankAndLocalRow[1][2].resize(3); // Rank 1 has 8 rows: local row 2 has 3 entries
    colValuePairsPerRankAndLocalRow[1][2][0] = ColValPair_t({1,0.0012461});
    colValuePairsPerRankAndLocalRow[1][2][1] = ColValPair_t({2,0.000160258});
    colValuePairsPerRankAndLocalRow[1][2][2] = ColValPair_t({4,-1.});

    colValuePairsPerRankAndLocalRow[1][3].resize(2); // Rank 1 has 8 rows: local row 3 has 2 entries
    colValuePairsPerRankAndLocalRow[1][3][0] = ColValPair_t({3,0.2});
    colValuePairsPerRankAndLocalRow[1][3][1] = ColValPair_t({6,-0.2});

    colValuePairsPerRankAndLocalRow[1][4].resize(2); // Rank 1 has 8 rows: local row 4 has 2 entries
    colValuePairsPerRankAndLocalRow[1][4][0] = ColValPair_t({2,-1.});
    colValuePairsPerRankAndLocalRow[1][4][1] = ColValPair_t({6,1.});

    colValuePairsPerRankAndLocalRow[1][5].resize(4); // Rank 1 has 8 rows: local row 5 has 4 entries
    colValuePairsPerRankAndLocalRow[1][5][0] = ColValPair_t({0,-9.29367e-05});
    colValuePairsPerRankAndLocalRow[1][5][1] = ColValPair_t({1,8.77169e-05});
    colValuePairsPerRankAndLocalRow[1][5][2] = ColValPair_t({5,5.21976e-06});
    colValuePairsPerRankAndLocalRow[1][5][3] = ColValPair_t({7,-1.});

    colValuePairsPerRankAndLocalRow[1][6].resize(4); // Rank 1 has 8 rows: local row 6 has 4 entries
    colValuePairsPerRankAndLocalRow[1][6][0] = ColValPair_t({3,-0.2});
    colValuePairsPerRankAndLocalRow[1][6][1] = ColValPair_t({4,1.});
    colValuePairsPerRankAndLocalRow[1][6][2] = ColValPair_t({6,0.2});
    colValuePairsPerRankAndLocalRow[1][6][3] = ColValPair_t({7,1.});

    colValuePairsPerRankAndLocalRow[1][7].resize(2); // Rank 1 has 8 rows: local row 7 has 2 entries
    colValuePairsPerRankAndLocalRow[1][7][0] = ColValPair_t({5,-1.});
    colValuePairsPerRankAndLocalRow[1][7][1] = ColValPair_t({6,1.});
  }
  else {
    throw std::runtime_error( "In TestCase2::prepareDataToCreateMatrix(): unsupported numRanks = " + std::to_string(this->m_numRanks) );
  }

  globalDomainIndicesPerRank = globalRowIndicesPerRank;
  globalRangeIndicesPerRank  = globalRowIndicesPerRank;

#if 0
Rank 0
tlpw_fromEpetra_afterSolverMapFwd_LHS: globalLength = 17, localLength = 9, numVectors = 1
m_tpetraMultiVector(v,i,value) = 0,0,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,1,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,2,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,3,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,4,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,5,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,6,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,7,0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,8,0.000000000000e+00
tlpw_fromEpetra_afterSolverMapFwd_RHS: globalLength = 17, localLength = 9, numVectors = 1
m_tpetraMultiVector(v,i,value) = 0,0,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,1,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,2,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,3,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,4,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,5,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,6,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,7,-0.000000000000e+00
m_tpetraMultiVector(v,i,value) = 0,8,-1.000000000000e+00
  
Rank 1
tlpw_fromEpetra_afterSolverMapFwd_LHS: globalLength = 17, localLength = 8, numVectors = 1
m_tpetraMultiVector(v,i,value) = 0,0,0
m_tpetraMultiVector(v,i,value) = 0,1,0
m_tpetraMultiVector(v,i,value) = 0,2,0
m_tpetraMultiVector(v,i,value) = 0,3,0
m_tpetraMultiVector(v,i,value) = 0,4,0
m_tpetraMultiVector(v,i,value) = 0,5,0
m_tpetraMultiVector(v,i,value) = 0,6,0
m_tpetraMultiVector(v,i,value) = 0,7,0
tlpw_fromEpetra_afterSolverMapFwd_RHS: globalLength = 17, localLength = 8, numVectors = 1
m_tpetraMultiVector(v,i,value) = 0,0,-2.31539e-05
m_tpetraMultiVector(v,i,value) = 0,1,-0
m_tpetraMultiVector(v,i,value) = 0,2,-0.000382438
m_tpetraMultiVector(v,i,value) = 0,3,-0
m_tpetraMultiVector(v,i,value) = 0,4,-0
m_tpetraMultiVector(v,i,value) = 0,5,2.31539e-05
m_tpetraMultiVector(v,i,value) = 0,6,-0
m_tpetraMultiVector(v,i,value) = 0,7,-0

#endif

}

#endif // TPETRA_CORE_TEST_REINDEX_TRANSFORM_TEST_CASE_2_DEF_HPP
