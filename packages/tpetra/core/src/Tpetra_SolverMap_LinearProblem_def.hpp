// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP
#define TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP

/// \file Tpetra_SolverMap_LinearProblem_def.hpp
/// \brief Definition of the Tpetra::SolverMap_LinearProblem class
///
/// If you want to use Tpetra::SolverMap_LinearProblem, include
/// "Tpetra_SolverMap_LinearProblem.hpp", a file which CMake generates
/// and installs for you).
///
/// If you only want the declaration of Tpetra::SolverMap_LinearProblem,
/// include "Tpetra_SolverMap_LinearProblem_decl.hpp".

#include <Tpetra_SolverMap_LinearProblem_decl.hpp>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::SolverMap_LinearProblem()
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~SolverMap_LinearProblem()
{
  if ((this->newObj_.get() != nullptr             ) &&
      (this->newObj_.get() != this->origObj_.get())) {
    this->newObj_.reset();
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()( const OriginalType & orig )
{
  this->origObj_ = orig;

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldMatrix = dynamic_cast< Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>* >(orig->getMatrix().get());
  Teuchos::RCP< MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > OldRHS    = orig->getRHS();
  Teuchos::RCP< MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> > OldLHS    = orig->getLHS();
  Teuchos::RCP< CrsMatrix  <Scalar, LocalOrdinal, GlobalOrdinal, Node> > NewMatrix = solverMapCrsMatrixTrans_( Teuchos::rcp< CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(OldMatrix) );

  if (NewMatrix.get() == OldMatrix) {
    // Same matrix, so use same problem
    this->newObj_ = this->origObj_;
  }
  else {
    this->newObj_ = Teuchos::rcp< LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(
                      new LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>( NewMatrix
                                                                                  , OldLHS
                                                                                  , OldRHS
                                                                                  ));
  }
  
  return this->newObj_;
}

//
// Explicit instantiation macro
//
// Must be expanded from within the Tpetra namespace!
//

#define TPETRA_SOLVERMAPLINEARPROBLEM_INSTANT(SCALAR,LO,GO,NODE) \
  template class SolverMap_LinearProblem< SCALAR , LO , GO , NODE >;

} // namespace Tpetra

#endif // TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP
