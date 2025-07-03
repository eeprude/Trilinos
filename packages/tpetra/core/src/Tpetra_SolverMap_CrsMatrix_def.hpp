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
//: NewGraph_ (nullptr)
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

  //if (NewGraph_) {
  //  std::cout << "EEP at solvermap_cm_transform<>::destructor(): calling NewGraph_.reset()" << std::endl;
  //  delete this->NewGraph_;
  //}
  
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
typename SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::construct( OriginalType orig )
{
  using map_t = Map<LocalOrdinal, GlobalOrdinal, Node>;

  assert( !orig.IndicesAreGlobal() );

  this->origObj_ = orig;

  // Test if matrix has missing local columns in its col map
  Teuchos::RCP<const map_t> RowMap    = orig->getRowMap();
  Teuchos::RCP<const map_t> DomainMap = orig->getDomainMap();
  Teuchos::RCP<const map_t> ColMap    = orig->getColMap();

  Teuchos::RCP<const Teuchos::Comm<int>> Comm = RowMap->getComm();
  size_t domain_localNumCols = DomainMap->getLocalNumElements();

  size_t localAmountOfDifferences(0);
  for (size_t i(0); i < domain_localNumCols; ++i) {
    if (DomainMap->getGlobalElement(i) != ColMap->getGlobalElement(i)) {
      localAmountOfDifferences += 1;
      break;
    }
  }
  size_t globalAmountOfDifferences(0);
  Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, 1, &localAmountOfDifferences, &globalAmountOfDifferences);
  std::cout << "In SolverMap_CrsMatrix<>::construct(): globalAmountOfDifferences = " << globalAmountOfDifferences << std::endl;
  
  if (false) { // (globalAmountOfDifferences == 0) { // AquiEEP
    this->newObj_ = this->origObj_;
  }
  else {
    using cg_t = CrsGraph<LocalOrdinal, GlobalOrdinal, Node>;
    using cm_t = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

    size_t localNumRows = RowMap->getLocalNumElements();

    std::cout << "Agora 000" << std::endl;
    
    // Create ColMap with all local rows included
    std::vector<GlobalOrdinal> globalColIndices(domain_localNumCols);

    // Fill globalColIndices list with GIDs of all local columns 
    for (size_t i(0); i < domain_localNumCols; ++i) {
      globalColIndices[i] = DomainMap->getGlobalElement(i);
    }

    std::cout << "Agora 001" << std::endl;

    // Now append to globalColIndices any ghost column entries
    size_t current_localNumCols = ColMap->getLocalNumElements();
    for (size_t i(0); i < current_localNumCols; ++i) {
      if ( DomainMap->isNodeGlobalElement( ColMap->getGlobalElement(i) ) ) {
        globalColIndices.push_back( ColMap->getGlobalElement(i) );
      }
    }
    
    std::cout << "Agora 002" << std::endl;

    size_t new_localNumCols = globalColIndices.size();
    size_t new_globalNumCols;
    Teuchos::reduceAll(*Comm, Teuchos::REDUCE_SUM, 1, &new_localNumCols, &new_globalNumCols);

    std::cout << "Agora 003" << std::endl;

    // Create new column map
    //NewColMap_ = new Epetra_Map( NewNumGlobalCols, NewNumMyCols,&globalColIndices[0], DomainMap.IndexBase64(), Comm );
    NewColMap_ = Teuchos::rcp<map_t>( new map_t( new_globalNumCols         // global_size_t       numGlobalElements
                                               , &globalColIndices[0]      // global_ordinal_type indexList[]
                                               , new_localNumCols          // local_ordinal_type  indexListSize
                                               , DomainMap->getIndexBase() // global_ordinal_type indexBase
                                               , Comm
                                               ) );

    std::cout << "Agora 004" << std::endl;

    // New Graph
    std::vector<size_t> NumIndicesPerRow(localNumRows);
    for (size_t i(0); i < localNumRows; ++i) {
      NumIndicesPerRow[i] = orig->getNumEntriesInLocalRow(i);
    }
    Teuchos::ArrayView<const size_t> array_NumIndicesPerRow(NumIndicesPerRow.data(), localNumRows);
    NewGraph_ = Teuchos::rcp<cg_t>( new cg_t( /*Teuchos::Copy,*/ RowMap, NewColMap_, array_NumIndicesPerRow ) );

    std::cout << "Agora 005" << std::endl;

    size_t MaxNumEntries = orig->getGlobalMaxNumRowEntries();
    std::vector<GlobalOrdinal> destIndices( MaxNumEntries );
    for (size_t i(0); i < localNumRows; ++i) {
      GlobalOrdinal RowGID = RowMap->getGlobalElement(i);
      std::cout << "Agora 006, i = " << i << std::endl;
      typename cg_t::nonconst_global_inds_host_view_type sourceIndices("si",MaxNumEntries);
      size_t numSourceIndices(0);
      //orig.Graph().ExtractGlobalRowCopy( RowGID, MaxNumEntries, NumEntries, &Indices[0] );
      orig->getGraph()->getGlobalRowCopy( RowGID, sourceIndices, numSourceIndices );
      std::cout << "Agora 007, i = " << i << std::endl;
      size_t NumEntries( numSourceIndices );
      for (size_t j(0); j < NumEntries; ++j) {
        destIndices[j] = sourceIndices[j];
      }
      //NewGraph_->InsertGlobalIndices( RowGID, NumEntries, &Indices[0] );
      NewGraph_->insertGlobalIndices( RowGID, NumEntries, destIndices.data() );
      std::cout << "Agora 008, i = " << i << std::endl;
    }
    std::cout << "Agora 009" << std::endl;
    Teuchos::RCP<const map_t> RangeMap = orig->getRangeMap();
    NewGraph_->fillComplete(DomainMap,RangeMap);

    std::cout << "Agora 010" << std::endl;

    // Initial construction of matrix 
    Teuchos::RCP<cm_t> NewMatrix = Teuchos::rcp<cm_t>( new cm_t( NewGraph_ ) );

    std::cout << "Agora 011" << std::endl;

#if 0
    // Insert views of row values
    int * myIndices;
    double * myValues;
    int indicesCnt;
#endif
    typename cm_t::local_inds_host_view_type orig_localRowIndices;
    typename cm_t::values_host_view_type     orig_localRowValues;

    typename cg_t::local_inds_host_view_type newGraph_localRowIndices;

    std::vector<Scalar>       newMatrix_localRowValues (MaxNumEntries);
    std::vector<LocalOrdinal> newMatrix_localRowIndices(MaxNumEntries);
    size_t new_localNumRows = NewMatrix->getLocalNumRows();
    std::cout << "Agora 012" << std::endl;
    for (size_t i(0); i < new_localNumRows; ++i) {
      //orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      orig->getLocalRowView( i
                           , orig_localRowIndices
                           , orig_localRowValues
                           );

      std::cout << "Agora 013, i = " << i << std::endl;

      //NewGraph_->ExtractMyRowView( i, indicesCnt, myIndices );
      NewGraph_->getLocalRowView( i
                                , newGraph_localRowIndices
                                );

      std::cout << "Agora 014, i = " << i << std::endl;

      assert( orig_localRowIndices.size() == newGraph_localRowIndices.size() );

      std::cout << "Agora 015, i = " << i << std::endl;

      //NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
      size_t numEntries( newGraph_localRowIndices.size() );
      for (size_t j(0); j < numEntries; ++j) {
        newMatrix_localRowValues [j] = orig_localRowValues [j];
        newMatrix_localRowIndices[j] = orig_localRowIndices[j];
      }
      NewMatrix->replaceLocalValues( i
                                   , numEntries
                                   , newMatrix_localRowValues.data()
                                   , newMatrix_localRowIndices.data()
                                   );
      std::cout << "Agora 016, i = " << i << std::endl;
    }

    std::cout << "Agora 017" << std::endl;

    NewMatrix->fillComplete(DomainMap,RangeMap);

    std::cout << "Agora 018" << std::endl;

    this->newObj_ = NewMatrix;
  }

  return this->newObj_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( const OriginalType & orig )
{
  return construct(orig);
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
