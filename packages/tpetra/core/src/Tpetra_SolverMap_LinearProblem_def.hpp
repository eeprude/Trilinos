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
/// \brief Definition of the Tpetra::LinearProblem_SolverMap class
///
/// If you want to use Tpetra::LinearProblem_SolverMap, include
/// "Tpetra_SolverMap_LinearProblem.hpp", a file which CMake generates
/// and installs for you).
///
/// If you only want the declaration of Tpetra::LinearProblem_SolverMap,
/// include "Tpetra_SolverMap_LinearProblem_decl.hpp".

#include <Tpetra_LinearProblem.h>
#include <Tpetra_CrsMatrix.h>

namespace Tpetra {

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
LinearProblem_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>::LinearProblem_SolverMap()
{
  // Nothing to do
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
LinearProblem_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>::~LinearProblem_SolverMap()
{
  if ((newObj_ != nullptr ) &&
      (newObj_ != origObj_)) {
    delete newObj_;
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename LinearProblem_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewTypeRef
LinearProblem_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()( typename LinearProblem_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node>::OriginalTypeRef orig )
{
  origObj_ = &orig;

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldMatrix = dynamic_cast< CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>* >( orig.GetMatrix() );
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldRHS = orig.GetRHS();
  MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldLHS = orig.GetLHS();

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> & NewMatrix = crsMatrixSolverMapTrans_( *OldMatrix );

  if( &NewMatrix == OldMatrix ) //same matrix so use same problem
    newObj_ = origObj_;
  else
    newObj_ = new LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>( &NewMatrix, OldLHS, OldRHS );

  return *newObj_;
}

} //namespace Tpetra

#endif // TPETRA_SOLVERMAP_LINEARPROBLEM_DEF_HPP
