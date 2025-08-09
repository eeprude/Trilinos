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
#include <Tpetra_Vector_decl.hpp>
#include <Tpetra_Import_decl.hpp>

#include <vector>

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
void
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::printMatrix( OriginalType const & passedMatrix ) const
{
  using cm_t = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  std::cout << ": globalNumRows = " << passedMatrix->getGlobalNumRows()
            << ", localNumRows = "  << passedMatrix->getLocalNumRows()
            << ", globalNumCols = " << passedMatrix->getGlobalNumCols()
            << ", localNumCols = "  << passedMatrix->getLocalNumCols()
            << ", globalNnz = "     << passedMatrix->getGlobalNumEntries()
            << ", localNnz = "      << passedMatrix->getLocalNumEntries()
            << ", tpetraMaxNnz = "  << passedMatrix->getGlobalMaxNumRowEntries()
            << std::endl;

  {
    std::cout << ": &rowMap() = " << passedMatrix->getRowMap().get()
              << ", rowMap()->globalNumElems = " << passedMatrix->getRowMap()->getGlobalNumElements()
              << ", localNumElems = " << passedMatrix->getRowMap()->getLocalNumElements()
              << "; globalIndices =";
    auto/*typename Map_t::global_indices_array_type*/ tmpArray = passedMatrix->getRowMap()->getMyGlobalIndices();
    for (size_t i(0); i < tmpArray.size(); ++i) {
      std::cout << " " << tmpArray[i];
    }
    std::cout << "; ";
    for (size_t i(0); i < passedMatrix->getRowMap()->getLocalNumElements(); ++i) {
      std::cout << " " << passedMatrix->getRowMap()->getGlobalElement(i);
    }
    std::cout << std::endl;
  }
  {
    std::cout << ": &colMap() = " << passedMatrix->getColMap().get()
              << ", colMap()->globalNumElems = " << passedMatrix->getColMap()->getGlobalNumElements()
              << ", localNumElems = " << passedMatrix->getColMap()->getLocalNumElements()
              << "; globalIndices =";
    auto/*typename Map_t::global_indices_array_type*/ tmpArray = passedMatrix->getColMap()->getMyGlobalIndices();
    for (size_t i(0); i < tmpArray.size(); ++i) {
      std::cout << " " << tmpArray[i];
    }
    std::cout << "; ";
    for (size_t i(0); i < passedMatrix->getColMap()->getLocalNumElements(); ++i) {
      std::cout << " " << passedMatrix->getColMap()->getGlobalElement(i);
    }
    std::cout << std::endl;
  }
  {
    std::cout << ": &domainMap() = " << passedMatrix->getDomainMap().get()
              << ", domainMap()->globalNumElems = " << passedMatrix->getDomainMap()->getGlobalNumElements()
              << ", localNumElems = " << passedMatrix->getDomainMap()->getLocalNumElements()
              << "; globalIndices =";
    auto/*typename Map_t::global_indices_array_type*/ tmpArray = passedMatrix->getDomainMap()->getMyGlobalIndices();
    for (size_t i(0); i < tmpArray.size(); ++i) {
      std::cout << " " << tmpArray[i];
    }
    std::cout << "; ";
    for (size_t i(0); i < passedMatrix->getDomainMap()->getLocalNumElements(); ++i) {
      std::cout << " " << passedMatrix->getDomainMap()->getGlobalElement(i);
    }
    std::cout << std::endl;
  }
  {
    std::cout << ": &rangeMap() = " << passedMatrix->getRangeMap().get()
              << ", rangeMap()->globalNumElems = " << passedMatrix->getRangeMap()->getGlobalNumElements()
              << ", localNumElems = " << passedMatrix->getRangeMap()->getLocalNumElements()
              << "; globalIndices =";
    auto/*typename Map_t::global_indices_array_type*/ tmpArray = passedMatrix->getRangeMap()->getMyGlobalIndices();
    for (size_t i(0); i < tmpArray.size(); ++i) {
      std::cout << " " << tmpArray[i];
    }
    std::cout << "; ";
    for (size_t i(0); i < passedMatrix->getRangeMap()->getLocalNumElements(); ++i) {
      std::cout << " " << passedMatrix->getRangeMap()->getGlobalElement(i);
    }
    std::cout << std::endl;
  }

  typename cm_t::values_host_view_type matrixValues;
  typename cm_t::local_inds_host_view_type matrixIndices;

  //int tpetraMaxNnz = passedMatrix->getGlobalMaxNumRowEntries();
  for (size_t i(0); i < passedMatrix->getLocalNumRows(); ++i) {
    int tpetraNnz = passedMatrix->getNumEntriesInLocalRow(i);
    assert( tpetraNnz <= tpetraMaxNnz);

    int tpetraNumEntries(tpetraNnz);
    passedMatrix->getLocalRowView( i
                                      , matrixIndices
                                      , matrixValues
                                      );

    std::cout << "i = " << i
              << ", tpetraNnz = " << tpetraNnz
              << ", numEntries = " << tpetraNumEntries
              << ";";
    std::cout << " indices =";
    for (size_t j(0); j < tpetraNumEntries; ++j) {
      std::cout << " " << matrixIndices[j];
    }
    std::cout << ";";
    for (size_t j(0); j < tpetraNumEntries; ++j) {
      std::cout << " (" << i << "," << matrixIndices[j] << ")=" << matrixValues[j];
    }
    std::cout << std::endl;
  }

  for (size_t i(0); i < passedMatrix->getLocalNumRows(); ++i) {
    int tpetraNnz = passedMatrix->getNumEntriesInLocalRow(i);
    assert( tpetraNnz <= tpetraMaxNnz);

    int tpetraNumEntries(tpetraNnz);
    passedMatrix->getLocalRowView( i
                                      , matrixIndices
                                      , matrixValues
                                      );

    if (i <= 9) {
      std::cout << " ";
    }
    std::cout << i << " t[";
    size_t j(0);
    for (size_t col(0); col < passedMatrix->getLocalNumCols(); ++col) {
      std::string tmp(" ");
      if (j < tpetraNumEntries) {
        if (col == matrixIndices[j]) {
          tmp = "X";
          j += 1;
        }
      }
      std::cout << " " << tmp;
    }
    std::cout << "]" << std::endl;
  }
}

template <class Scalar,
          class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
typename Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::NewType
Reindex_CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>::operator()( OriginalType const & origMatrix )
{
  Teuchos::RCP< const Teuchos::Comm<int> > const & teuchosComm = origMatrix->getDomainMap()->getComm();
  int myRank( teuchosComm->getRank() );
  int numRanks( teuchosComm->getSize() );

  std::cout << "Entering Tpetra Reindex_CrsMatrix<>::operator()" << std::endl;
  using cm_t = CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>;

  //this->printMatrix(origMatrix);

  this->origObj_ = origMatrix;

  assert( origMatrix->getRowMap()->getLocalNumElements() == newRowMap_->getLocalNumElements() );

  // AquiEEP orig.OperatorDomainMap()
  if ((origMatrix->getDomainMap()->getGlobalNumElements() == 0) &&
      (origMatrix->getRowMap()->getGlobalNumElements()    == 0)) {
    std::cout << "In Tpetra Reindex_CrsMatrix<>::operator(): zero matrix situation" << std::endl;
    // Construct a zero matrix as a placeholder, don't do reindexing analysis.
    this->newObj_ = Teuchos::rcp<cm_t>( new cm_t(origMatrix->getRowMap(), origMatrix->getColMap(), 0) );
  }
  else {
    using map_t = Map   <               LocalOrdinal, GlobalOrdinal, Node>;
    using imp_t = Import<               LocalOrdinal, GlobalOrdinal, Node>;
    using v_t   = Vector<GlobalOrdinal, LocalOrdinal, GlobalOrdinal, Node>;

    // Construct new column map
    v_t cols( origMatrix->getDomainMap() );
    { 
      size_t origDomainMap_localSize = origMatrix->getDomainMap()->getLocalNumElements();
      map_t tmpColMap( origMatrix->getDomainMap()->getGlobalNumElements(), origDomainMap_localSize, 0, origMatrix->getDomainMap()->getComm() );
      for (size_t i(0); i < origDomainMap_localSize; ++i) {
        cols.replaceLocalValue(i, tmpColMap.getGlobalElement(i));
      }
    }
    {
      std::cout.flush();
      teuchosComm->barrier();
      for (int p(0); p < numRanks; ++p) {
        teuchosComm->barrier();
        if (p != myRank) continue;

        using size_type = typename Teuchos::ArrayView<Scalar>::size_type;
        Teuchos::ArrayRCP<const GlobalOrdinal> aux = cols.getData();
        std::cout << "In Tpetra Reindex_CrsMatrix<>::operator()"
                  << ": cols =";
        for (size_type j(0); j < aux.size(); ++j) {
          std::cout << " " << aux[j];
        }
        std::cout << std::endl;

	std::cout.flush();
      }
      teuchosComm->barrier();
    }

    imp_t importer( origMatrix->getDomainMap(), origMatrix->getColMap() );
    v_t newCols( origMatrix->getColMap() );
    newCols.doImport( cols     // const SrcDistObject& source
                    , importer // const Import<LocalOrdinal, GlobalOrdinal, Node>& importer
                    , INSERT   // const CombineMode CM
                    , false    // const bool restrictedMode
                    );
    {
      using size_type = typename Teuchos::ArrayView<Scalar>::size_type;
      Teuchos::ArrayRCP<const GlobalOrdinal> aux = newCols.getData();
      std::cout << "In Tpetra Reindex_CrsMatrix<>::operator()"
                << ": newCols =";
      for (size_type j(0); j < aux.size(); ++j) {
        std::cout << " " << aux[j];
      }
      std::cout << std::endl;
    }

    Teuchos::ArrayRCP<const GlobalOrdinal> newColIndicesArray = newCols.getData();
    std::vector<GlobalOrdinal> newColIndicesVector(newColIndicesArray.size());
    for (size_t j(0); j < newColIndicesVector.size(); ++j) {
      newColIndicesVector[j] = newColIndicesArray[j];
    }
    this->newColMap_ = Teuchos::RCP<map_t>( new map_t( origMatrix->getColMap()->getGlobalNumElements() // const global_size_t numGlobalElements
                                                     , newColIndicesVector.data()                      // const global_ordinal_type indexList[]
                                                     , origMatrix->getColMap()->getLocalNumElements()  // const local_ordinal_type indexListSize
                                                     , origMatrix->getColMap()->getIndexBase()         // const global_ordinal_type indexBase
                                                     , origMatrix->getColMap()->getComm()              // const Teuchos::RCP<const Teuchos::Comm<int> >& comm
                                                     ));

    {
      auto aux = this->newRowMap_->getMyGlobalIndices();
      std::cout << "In Tpetra Reindex_CrsMatrix<>::operator()"
                << ": newRowMap_ =";
      for (size_t j(0); j < aux.size(); ++j) {
        std::cout << " " << aux[j];
      }
      std::cout << std::endl;
    }
    {
      auto aux = this->newColMap_->getMyGlobalIndices();
      std::cout << "In Tpetra Reindex_CrsMatrix<>::operator()"
                << ": newColMap_ =";
      for (size_t j(0); j < aux.size(); ++j) {
        std::cout << " " << aux[j];
      }
      std::cout << std::endl;
    }
    //this->printMatrix(origMatrix);

    // Create the new matrix 
    size_t const origMatrix_maxNumEntries = origMatrix->getGlobalMaxNumRowEntries();
    Teuchos::RCP<cm_t> newMatrix = Teuchos::rcp<cm_t>( new cm_t(this->newRowMap_, this->newColMap_, origMatrix_maxNumEntries) );

    std::vector<Scalar>       newMatrix_localValues (origMatrix_maxNumEntries);
    std::vector<LocalOrdinal> newMatrix_localIndices(origMatrix_maxNumEntries);

    typename cm_t::local_inds_host_view_type origMatrix_localIndices;
    typename cm_t::values_host_view_type     origMatrix_localValues;

    size_t const newMatrix_localNumRows = newMatrix->getLocalNumRows();
    for (size_t i(0); i < newMatrix_localNumRows; ++i) {
      //orig.ExtractMyRowView( i, indicesCnt, myValues, myIndices );
      //NewMatrix->InsertMyValues( i, indicesCnt, myValues, myIndices );

      origMatrix->getLocalRowView( i, origMatrix_localIndices, origMatrix_localValues );

      size_t const numEntries( origMatrix_localIndices.size() );
      for (size_t j(0); j < numEntries; ++j) {
        newMatrix_localValues [j] = origMatrix_localValues[j];
        newMatrix_localIndices[j] = origMatrix_localIndices[j];
      }

      std::cout << "In Tpetra Reindex_CrsMatrix<>::operator()"
                << ", calling newMatrix->replaceLocalValues()"
                << ": i = " << i
                << ", numEntries = " << numEntries
                << ", newMatrix_localValues =";
      for (size_t j(0); j < numEntries; ++j) {
        std::cout << " " << newMatrix_localValues[j];
      }
      std::cout << ", newMatrix_localIndices =";
      for (size_t j(0); j < numEntries; ++j) {
        std::cout << " " << newMatrix_localIndices[j];
      }
      std::cout << std::endl;
      newMatrix->insertLocalValues( i                             // const LocalOrdinal localRow
                                  , numEntries                    // const LocalOrdinal numEnt
                                  , newMatrix_localValues.data()  // const Scalar       vals[]
                                  , newMatrix_localIndices.data() // const LocalOrdinal cols[]
                                  , INSERT                        // const CombineMode  CM
                                  );
    }

    newMatrix->fillComplete();

    this->newObj_ = newMatrix;
  }

  std::cout << "Leaving Tpetra Reindex_CrsMatrix<>::operator()" << std::endl;

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

