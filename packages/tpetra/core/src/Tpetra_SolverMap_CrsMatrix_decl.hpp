// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef TPETRA_SOLVERMAP_CRSMATRIX_DECL_HPP
#define TPETRA_SOLVERMAP_CRSMATRIX_DECL_HPP

/// \file Tpetra_SolverMap_CrsMatrix_decl.hpp
/// \brief Declaration of the Tpetra::CrsMatrix_SolverMap class

#include <Tpetra_Transform.h>
#include <Tpetra_CrsMatrix.h>

namespace Tpetra {

///
/** Given an input CrsMatrix, the column map is checked for missing indices associated
 *  with the local rows.  If found, a view of the CrsMatrix is formed using a new
 *  CrsGraph with a fixed column mapping including all local row indices.
 */
template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class CrsMatrix_SolverMap<Scalar, LocalOrdinal, GlobalOrdinal, Node> : public StructuralSameTypeTransform< CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
{
public:
  ///
  /** Constructor
   */
  CrsMatrix_SolverMap();

  ///
  /** Destructor
   */
  ~CrsMatrix_SolverMap();

  ///
  /** Constructs fixed view of CrsMatrix as necessary.
   */
  NewTypeRef operator()( OriginalTypeRef orig );

private:
  template<typename int_type>
  NewTypeRef construct( OriginalTypeRef orig );

  // avoid virtual function hidden warning
  using StructuralSameTypeTransform< CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >::construct;

  Map<LocalOrdinal, GlobalOrdinal, Node> * NewColMap_;
  CrsGraph<LocalOrdinal, GlobalOrdinal, Node> * NewGraph_;
};

} // namespace Tpetra

#endif // TPETRA_SOLVERMAP_CRSMATRIX_DECL_HPP
