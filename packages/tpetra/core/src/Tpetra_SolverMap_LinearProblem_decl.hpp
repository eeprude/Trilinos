// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_LINEARPROBLEM_DECL_HPP
#define TPETRA_SOLVERMAP_LINEARPROBLEM_DECL_HPP

/// \file Tpetra_SolverMap_LinearProblem_decl.hpp
/// \brief Declaration of the Tpetra::LinearProblem_SolverMap class

#include <Tpetra_Transform.h>

#include <Tpetra_SolverMap_CrsMatrix.h>

namespace Tpetra {

///
/** Constructs a LinearProblem with a "fixed" Column Map for the CrsMatrix.
 *  Almost entirely a view except for the "fixed" Tpetra_CrsGraph.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class LinearProblem_SolverMap : public StructuralSameTypeTransform< LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
{
public:
  ///
  /** Constructor
   */
  LinearProblem_SolverMap();

  ///
  /** Destructor
   */
  ~LinearProblem_SolverMap();

  ///
  /** Constructs "fixed" Tpetra::LinearProblem
   */
  NewTypeRef operator()( OriginalTypeRef orig );

private:
  CrsMatrix_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node> crsMatrixSolverMapTrans_;
};

} //namespace Tpetra

#endif // TPETRA_SOLVERMAP_LINEARPROBLEM_DECL_HPP
