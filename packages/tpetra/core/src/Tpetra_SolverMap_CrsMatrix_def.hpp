// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_CRSMATRIX_DEF_HPP
#define TPETRA_SOLVERMAP_CRSMATRIX_DEF_HPP

/// \file Tpetra_SolverMap_CrsMatrix_def.hpp
/// \brief Definition of the Tpetra::SolverMap_CrsMatrix class
///
/// If you want to use Tpetra::SolverMap_CrsMatrix, include
/// "Tpetra_SolverMap_CrsMatrix.hpp", a file which CMake generates
/// and installs for you).
///
/// If you only want the declaration of Tpetra::SolverMap_CrsMatrix,
/// include "Tpetra_SolverMap_CrsMatrix_decl.hpp".

#include <Tpetra_SolverMap_CrsMatrix_decl.hpp>

#include <vector>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SolverMap_CrsMatrix()
//: NewColMap_(nullptr)
  : NewGraph_ (nullptr)
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SolverMap_CrsMatrix()
{
  std::cout << "EEP at solvermap_cm_transform<>::destructor(): entering" << std::endl;
  if ((this->newObj_.get() != nullptr             ) &&
      (this->newObj_.get() != this->origObj_.get())) {
    std::cout << "EEP at solvermap_cm_transform<>::destructor(): calling newObj_.reset()" << std::endl;
    this->newObj_.reset();
  }

  if (NewGraph_) {
    std::cout << "EEP at solvermap_cm_transform<>::destructor(): calling NewGraph_.reset()" << std::endl;
    delete this->NewGraph_;
  }
  
  //if (NewColMap_) {
  //  std::cout << "EEP at solvermap_cm_transform<>::destructor(): calling NewColMap_.reset()" << std::endl;
  //  delete this->NewColMap_;
  //}
  std::cout << "EEP at solvermap_cm_transform<>::destructor(): leaving" << std::endl;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
template<typename int_type>
typename SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::construct( OriginalType orig )
{
  using map_t = Map<LocalOrdinal, GlobalOrdinal, Node>;
  using cg_t = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;

  this->origObj_ = orig;

  assert( !orig.IndicesAreGlobal() );
#if 0 // AquiEEP
  this->newObj_ = this->origObj_;
#else
  //test if matrix has missing local columns in its col std::map
  Teuchos::RCP<const map_t> RowMap    = orig->getRowMap();
  Teuchos::RCP<const map_t> DomainMap = orig->getDomainMap();
  Teuchos::RCP<const map_t> ColMap    = orig->getColMap();

  Teuchos::RCP<const Teuchos::Comm<int>> Comm = RowMap->getComm();
  int localNumRows = RowMap->getLocalNumElements();
  size_t domain_localNumCols = DomainMap->getLocalNumElements();

  size_t localNumDifferences = 0;
  for (size_t i(0); i < domain_localNumCols; ++i) {
    if (DomainMap->getGlobalElement(i) != ColMap->getGlobalElement(i)) {
      localNumDifferences = 1;
      break;
    }
  }

  size_t globalNumDifferences = 0;
  Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, 1, &localNumDifferences, &globalNumDifferences);
  
  if (globalNumDifferences == 0) {
    this->newObj_ = this->origObj_;
  }
  else {
    // Create ColMap with all local rows included
    std::vector<GlobalOrdinal> Cols(domain_localNumCols);

    // Fill Cols list with GIDs of all local columns 
    for (size_t i(0); i < domain_localNumCols; ++i) {
      Cols[i] = DomainMap->getGlobalElement(i);
    }

    // Now append to Cols any ghost column entries
    int current_localNumCols = ColMap->getLocalNumElements();
    for (int i(0); i < current_localNumCols; ++i) {
      if ( DomainMap->isNodeGlobalElement( ColMap->getGlobalElement(i) ) ) {
        Cols.push_back( ColMap->getGlobalElement(i) );
      }
    }
    
    size_t new_localNumCols = Cols.size();
    size_t new_globalNumCols;
    Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, 1, &new_localNumCols, &new_globalNumCols);

    // Create new column std::map
    NewColMap_ = Teuchos::rcp<map_t>( new map_t( new_globalNumCols, new_localNumCols, /*&Cols[0],*/ DomainMap->getIndexBase(), Comm ) ); // AquiEEP

    // New Graph
    std::vector<size_t> NumIndicesPerRow(localNumRows);
    for (int i(0); i < localNumRows; ++i) {
      NumIndicesPerRow[i] = orig->getNumEntriesInLocalRow(i);
    }
    Teuchos::ArrayView<const size_t> array_NumIndicesPerRow(NumIndicesPerRow.data(), localNumRows);
    NewGraph_ = new cg_t( /*Teuchos::Copy,*/ RowMap, NewColMap_, array_NumIndicesPerRow );

    size_t MaxNumEntries = orig->getGlobalMaxNumRowEntries();
    //size_t NumEntries;
    std::vector<GlobalOrdinal> Indices( MaxNumEntries );
    for (int i(0); i < localNumRows; ++i) {
      //GlobalOrdinal RowGID = RowMap->getGlobalElement(i);
      //orig.Graph().ExtractGlobalRowCopy( RowGID, MaxNumEntries, NumEntries, &Indices[0] );
      //NewGraph_->InsertGlobalIndices( RowGID, NumEntries, &Indices[0] );
    }
#if 0
    const Epetra_Map & RangeMap = orig.RangeMap();
    NewGraph_->FillComplete(DomainMap,RangeMap);

    // Intial construction of matrix 
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, *NewGraph_ );

    // Insert views of row values
    int * myIndices;
    double * myValues;
    int indicesCnt;
    int numMyRows = NewMatrix->localNumRows();
    for( int i = 0; i < numMyRows; ++i )
    {
      orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      NewGraph_->ExtractMyRowView( i, indicesCnt, myIndices );

      NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
    }

    NewMatrix->FillComplete(DomainMap,RangeMap);

    this->newObj_ = NewMatrix;
#endif
  }
#endif
  return this->newObj_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( const OriginalType & orig )
{
  return construct<GlobalOrdinal>(orig); // AquiEEP
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_SOLVERMAPCRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  template class SolverMap_CrsMatrix< SCALAR , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_SOLVERMAP_CRSMATRIX_DEF_HPP
