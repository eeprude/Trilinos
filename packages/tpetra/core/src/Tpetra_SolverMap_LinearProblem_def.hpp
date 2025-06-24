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
  if ((this->newObj_ != nullptr       ) &&
      (this->newObj_ != this->origObj_)) {
    delete this->newObj_;
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewTypeRef
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()( typename SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::OriginalTypeRef orig )
{
  this->origObj_ = &orig;

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldMatrix = dynamic_cast< CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>* >( orig.GetMatrix() );
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldRHS = orig.GetRHS();
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldLHS = orig.GetLHS();

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & NewMatrix = solverMapCrsMatrixTrans_( *OldMatrix );

  if( &NewMatrix == OldMatrix ) //same matrix so use same problem
    this->newObj_ = this->origObj_;
  else
    this->newObj_ = new LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>( &NewMatrix, OldLHS, OldRHS );

  return *(this->newObj_);
}

} // namespace Tpetra

#endif // TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP
