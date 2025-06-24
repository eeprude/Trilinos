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
  : NewColMap_(nullptr)
  , NewGraph_ (nullptr)
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SolverMap_CrsMatrix()
{
  if ((this->newObj_.get() != nullptr             ) &&
      (this->newObj_.get() != this->origObj_.get())) {
    this->newObj_.reset();
  }

  if (NewGraph_) {
    delete this->NewGraph_;
  }
  
  if (NewColMap_) {
    delete this->NewColMap_;
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
template<typename int_type>
typename SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::construct( OriginalType orig )
{
  this->origObj_ = orig;

  assert( !orig.IndicesAreGlobal() );
#if 0 // AquiEEP
  //test if matrix has missing local columns in its col std::map
  const Epetra_Map & RowMap = orig.RowMap();
  const Epetra_Map & DomainMap = orig.DomainMap();
  const Epetra_Map & ColMap = orig.ColMap();
  const Epetra_Comm & Comm = RowMap.Comm();
  int NumMyRows = RowMap.NumMyElements();
  int NumCols = DomainMap.NumMyElements();
  int Match = 0;
  for( int i = 0; i < NumCols; ++i )
    if( DomainMap.GID64(i) != ColMap.GID64(i) )
    {
      Match = 1;
      break;
    }

  int MatchAll = 0;
  Comm.SumAll( &Match, &MatchAll, 1 );

  if( !MatchAll )
  {
    this->newObj_ = this->origObj_;
  }
  else
  {
    //create ColMap with all local rows included
    std::vector<int_type> Cols(NumCols);
    //fill Cols list with GIDs of all local columns 
    for( int i = 0; i < NumCols; ++i )
      Cols[i] = (int_type) DomainMap.GID64(i);

    //now append to Cols any ghost column entries
    int NumMyCols = ColMap.NumMyElements();
    for( int i = 0; i < NumMyCols; ++i )
      if( !DomainMap.MyGID( ColMap.GID64(i) ) ) Cols.push_back( (int_type) ColMap.GID64(i) );
    
    int NewNumMyCols = Cols.size();
    int NewNumGlobalCols;
    Comm.SumAll( &NewNumMyCols, &NewNumGlobalCols, 1 );
    //create new column std::map
    NewColMap_ = new Epetra_Map( NewNumGlobalCols, NewNumMyCols,&Cols[0], DomainMap.IndexBase64(), Comm );

    //New Graph
    std::vector<int> NumIndicesPerRow( NumMyRows );
    for( int i = 0; i < NumMyRows; ++i )
      NumIndicesPerRow[i] = orig.NumMyEntries(i);
    NewGraph_ = new Epetra_CrsGraph( Copy, RowMap, *NewColMap_, &NumIndicesPerRow[0] );

    int MaxNumEntries = orig.MaxNumEntries();
    int NumEntries;
    std::vector<int_type> Indices( MaxNumEntries );
    for( int i = 0; i < NumMyRows; ++i )
    {
      int_type RowGID = (int_type) RowMap.GID64(i);
      orig.Graph().ExtractGlobalRowCopy( RowGID, MaxNumEntries, NumEntries, &Indices[0] );
      NewGraph_->InsertGlobalIndices( RowGID, NumEntries, &Indices[0] );
    }
    const Epetra_Map & RangeMap = orig.RangeMap();
    NewGraph_->FillComplete(DomainMap,RangeMap);

    //intial construction of matrix 
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, *NewGraph_ );

    //insert views of row values
    int * myIndices;
    double * myValues;
    int indicesCnt;
    int numMyRows = NewMatrix->NumMyRows();
    for( int i = 0; i < numMyRows; ++i )
    {
      orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      NewGraph_->ExtractMyRowView( i, indicesCnt, myIndices );

      NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
    }

    NewMatrix->FillComplete(DomainMap,RangeMap);

    this->newObj_ = NewMatrix;
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
#if 0 // AquiEEP
#ifndef EPETRA_NO_32BIT_GLOBAL_INDICES
  if(orig.RowMap().GlobalIndicesInt()) {
    return construct<int>(orig);
  }
  else
#endif
#ifndef EPETRA_NO_64BIT_GLOBAL_INDICES
  if(orig.RowMap().GlobalIndicesLongLong()) {
    return construct<long long>(orig);
  }
  else
#endif
    throw "SolverMap_CrsMatrix::operator(): GlobalIndices type unknown";
#endif
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
