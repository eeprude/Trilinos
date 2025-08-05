// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_REINDEX_MULTIVECTOR_DEF_HPP
#define TPETRA_REINDEX_MULTIVECTOR_DEF_HPP

#include <Tpetra_Reindex_MultiVector_decl.hpp>

/// \file Tpetra_Reindex_MultiVector_def.hpp
/// \brief Definition of the Tpetra::Reindex_MultiVector class
///
/// If you want to use Tpetra::Reindex_MultiVector, include
/// "Tpetra_Reindex_MultiVector.hpp", a file which CMake generates
/// and installs for you.
///
/// If you only want the declaration of Tpetra::Reindex_MultiVector,
/// include "Tpetra_Reindex_MultiVector_decl.hpp".

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
Reindex_MultiVector( Teuchos::RCP< Map<LocalOrdinal, GlobalOrdinal, Node> const > newRowMap )
  : ViewTransform< MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> >()
  , newRowMap_(newRowMap)
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~Reindex_MultiVector()
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( OriginalType const & orig )
{
  this->origObj_ = orig;
#if 0 // Aqui
  //test std::map, must have same number of local and global elements as original row std::map
  assert( orig.Map().NumMyElements() == newRowMap_.NumMyElements() );

  std::vector<double*> MyValues(1);
  int MyLDA;
  int NumVectors = orig.NumVectors();
  orig.ExtractView( &MyValues[0], &MyLDA );

  Epetra_MultiVector * newMV = new Epetra_MultiVector( View, newRowMap_, MyValues[0], MyLDA, NumVectors );

  this->newObj_ = newMV;
#endif
  return this->newObj_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::transform( OriginalType orig )
{
  //test std::map, must have same number of local and global elements as original row std::map
  assert( orig.Map().NumMyElements() == newRowMap_.NumMyElements() );
#if 0 // Aqui
  std::vector<double*> MyValues(1);
  int MyLDA;
  int NumVectors = orig.NumVectors();
  orig.ExtractView( &MyValues[0], &MyLDA );

  return Teuchos::rcp(new Epetra_MultiVector( View, newRowMap_, MyValues[0], MyLDA, NumVectors ));
#else 
  return this->newObj_;
#endif
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_REINDEXMULTIVECTOR_INSTANT(SCALAR,LO,GO,NODE) \
  template class Reindex_MultiVector< SCALAR , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_REINDEX_MULTIVECTOR_DEF_HPP
