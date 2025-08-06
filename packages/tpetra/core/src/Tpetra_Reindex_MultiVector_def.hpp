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
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( OriginalType const & origMultiVector )
{
  using mv_t = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  this->origObj_ = origMultiVector;
  assert( origMultiVector->getMap()->getLocalNumElements() == this->newRowMap_->getLocalNumElements() );
  assert( origMultiVector->isConstantStride() == true );

  Teuchos::ArrayRCP<Teuchos::ArrayRCP<const Scalar> > origValues = origMultiVector->get2dView();

  using size_type = typename Teuchos::ArrayRCP<Scalar>::size_type;
  size_type numEntries(0);
  for (size_type j(0); j < origValues.size(); ++j) {
    numEntries += origValues[j].size();
  }

  std::vector<Scalar> tmpVec(numEntries);
  size_t k(0);
  for (size_type j(0); j < origValues.size(); ++j) {
    for (size_type i(0); i < origValues[j].size(); ++i) {
      tmpVec[k++] = origValues[j][i];
    }
  }

  Teuchos::ArrayView<Scalar const> valuesToInsert(tmpVec.data(), numEntries);
  
  this->newObj_ = Teuchos::RCP<mv_t>( new mv_t( newRowMap_                       // const Teuchos::RCP<const map_type>& map,
                                              , valuesToInsert                   // const Teuchos::ArrayView<const Scalar>& A,
                                              , origMultiVector->getStride()     // const size_t LDA,
                                              , origMultiVector->getNumVectors() // const size_t NumVectors
                                              ));

  return this->newObj_;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>::transform( OriginalType origMultiVector )
{
  using mv_t = MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  throw std::runtime_error("Reindex_MultiVector<>::transform() is not implemented yet"); // Aqui
  assert( origMultiVector->getMap()->getLocalNumElements() == this->newRowMap_->getLocalNumElements() );

#if 0 // Aqui
  std::vector<double*> MyValues(1);
  int MyLDA;
  int NumVectors = origMultiVector.NumVectors();
  origMultiVector.ExtractView( &MyValues[0], &MyLDA );

  return Teuchos::rcp(new mv_t( View, newRowMap_, MyValues[0], MyLDA, NumVectors ));
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
