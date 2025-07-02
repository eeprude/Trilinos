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
  std::cout << "EEP at solvermap_lp_transform<>::destructor(): entering" << std::endl;
  //if ((this->newObj_.get() != nullptr             ) &&
  //    (this->newObj_.get() != this->origObj_.get())) {
  //  std::cout << "EEP at solvermap_lp_transform<>::destructor(): calling newObj_.reset()" << std::endl;
  //  this->newObj_.reset();
  //}
  std::cout << "EEP at solvermap_lp_transform<>::destructor(): leaving" << std::endl;
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
SolverMap_LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
operator()( const OriginalType & orig )
{
  using mv_t = MultiVector  <Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using cm_t = CrsMatrix    <Scalar, LocalOrdinal, GlobalOrdinal, Node>;
  using lp_t = LinearProblem<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  std::cout << "EEP at solverMap_lp_transform::operator(): entering" << std::endl;

  this->origObj_ = orig;

  CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> * OldMatrix = dynamic_cast< Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>* >(orig->getMatrix().get());
  Teuchos::RCP< mv_t > OldRHS    = orig->getRHS();
  Teuchos::RCP< mv_t > OldLHS    = orig->getLHS();
  std::cout << "EEP at solverMap_lp_transform::operator(): calling solverMapCrsMatrixTrans_()" << std::endl;
  Teuchos::RCP< cm_t > NewMatrix = solverMapCrsMatrixTrans_( Teuchos::rcp< cm_t >(OldMatrix) );
  std::cout << "EEP at solverMap_lp_transform::operator(): returned from solverMapCrsMatrixTrans_()" << std::endl;

  if (false) { // (NewMatrix.get() == OldMatrix) { // AquiEEP
    // Same matrix, so use same problem
    std::cout << "EEP at solverMap_lp_transform::operator(): simple =" << std::endl;
    this->newObj_ = this->origObj_;
  }
  else {
    std::cout << "EEP at solverMap_lp_transform::operator(): new lp_t" << std::endl;
    this->newObj_ = Teuchos::rcp< lp_t >( new lp_t(NewMatrix, OldLHS, OldRHS) );
  }
  
  std::cout << "EEP at solverMap_lp_transform::operator(): leaving" << std::endl;
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
