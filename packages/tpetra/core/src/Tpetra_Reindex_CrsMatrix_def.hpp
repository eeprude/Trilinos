// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_CRSMATRIX_DEF_HPP
#define TPETRA_REINDEX_CRSMATRIX_DEF_HPP

#include <Tpetra_Reindex_CrsMatrix_decl.hpp>

#include <vector>

//#include <Epetra_Export.h> // Aqui
//#include <Epetra_Import.h>

/// \file Tpetra_Reindex_CrsMatrix_def.hpp
/// \brief Definition of the Tpetra::Reindex_CrsMatrix class
///
/// If you want to use Tpetra::Reindex_CrsMatrix, include
/// "Tpetra_Reindex_CrsMatrix.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::Reindex_CrsMatrix,
/// include "Tpetra_Reindex_CrsMatrix_decl.hpp".

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Reindex_CrsMatrix( Teuchos::RCP< Map<LocalOrdinal, GlobalOrdinal, Node> const > newRowMap )
  : ViewTransform< CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , newRowMap_(newRowMap)
  , newColMap_(Teuchos::null)
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Reindex_CrsMatrix()
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( OriginalType const & orig )
{
  this->origObj_ = orig;
#if 0 // Aqui
  //test std::map, must have same number of local and global elements as original row std::map
  //Epetra_Map & OldRowMap = const_cast<Epetra_Map&>(orig.RowMap());
  Epetra_Map & OldDomainMap = const_cast<Epetra_Map&>(orig.OperatorDomainMap());
  Epetra_Map & OldColMap = const_cast<Epetra_Map&>(orig.ColMap());
  int NumMyElements = OldDomainMap.NumMyElements();
  int_type NumGlobalElements = (int_type) OldDomainMap.NumGlobalElements64();
  assert( orig.RowMap().NumMyElements() == NewRowMap_.NumMyElements() );

  if (NumGlobalElements == 0 && orig.RowMap().NumGlobalElements64() == 0 )
  {
    //construct a zero matrix as a placeholder, don't do reindexing analysis.
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, orig.RowMap(), orig.ColMap(), 0 );
    newObj_ = NewMatrix;
  }
  else {

    //Construct new Column Map
    typename Epetra_GIDTypeVector<int_type>::impl Cols( OldDomainMap );
    typename Epetra_GIDTypeVector<int_type>::impl NewCols( OldColMap );
    Epetra_Import Importer( OldColMap, OldDomainMap );
 
    Epetra_Map tmpColMap( NumGlobalElements, NumMyElements, 0, OldDomainMap.Comm() );
 
    for( int i = 0; i < NumMyElements; ++i )
      Cols[i] = (int_type) tmpColMap.GID64(i);

    NewCols.Import( Cols, Importer, Insert );

    std::vector<int_type*> NewColIndices(1);
    NewCols.ExtractView( &NewColIndices[0] );

    int NumMyColElements = OldColMap.NumMyElements();
    int_type NumGlobalColElements = (int_type) OldColMap.NumGlobalElements64();

    NewColMap_ = new Epetra_Map( NumGlobalColElements, NumMyColElements, NewColIndices[0], (int_type) OldColMap.IndexBase64(), OldColMap.Comm() );

    //intial construction of matrix 
    Epetra_CrsMatrix * NewMatrix = new Epetra_CrsMatrix( View, NewRowMap_, *NewColMap_, 0 );

    //insert views of row values
    int * myIndices;
    double * myValues;
    int indicesCnt;
    int numMyRows = NewMatrix->NumMyRows();
    for( int i = 0; i < numMyRows; ++i )
    {
      orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );
    }

    NewMatrix->FillComplete();

    newObj_ = NewMatrix;
  }
#endif
  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REINDEXCRSMATRIX_INSTANT(SCALAR,LO,GO,NODE) \
  template class Reindex_CrsMatrix< SCALAR , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_REINDEX_CRSMATRIX_DEF_HPP

