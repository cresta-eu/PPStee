/*! \file mainSystemTest.cpp
 * \brief Main system test file.
 *
 * \date 14.05.2013
 * \author matu_gr
 *
 * \copyright Copyright 2013 German Aerospace Center (http://www.DLR.de)\n\n
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 */

#include "gtest/gtest.h"

#include "ppstee.hpp"

#include <unistd.h>

#define VERBOSE 1

#define COUTMPIME cout << "[" << mpiMe << "] "
#define OSMPIME os << "[" << pgraph->getMpiMe() << "] "

using std::cout;
using std::endl;
using std::ostream;

using namespace testing;

// include graph data
#ifdef USING_NUMPROCS_2
  #include "test/res/graphData-2.cpp"
#endif
#ifdef USING_NUMPROCS_3
  #include "test/res/graphData-3.cpp"
#endif
#ifdef USING_NUMPROCS_4
  #include "test/res/graphData-4.cpp"
#endif
#ifdef USING_NUMPROCS_5
  #include "test/res/graphData-5.cpp"
#endif
#ifdef USING_NUMPROCS_6
  #include "test/res/graphData-6.cpp"
#endif
#ifdef USING_NUMPROCS_12
  #include "test/res/graphData-12.cpp"
#endif
#ifdef USING_NUMPROCS_113
  #include "test/res/graphData-113.cpp"
#endif

class GraphDataProvider {
  protected:
    static int mpiN;
    static int mpiMe;

    static int* exVertglbtab;
    static int* exVertloctab;
    static int* exEdgeloctab;
    static int* exVertwgtcmp;
    static int* exEdgewgtcmp;
    static int* exVertwgtvis;
    static int* exEdgewgtvis;

    static void initGraphData() {
      MPI_Comm_size(MPI_COMM_WORLD, &GraphDataProvider::mpiN);
      MPI_Comm_rank(MPI_COMM_WORLD, &GraphDataProvider::mpiMe);

      // get example graph data
      exGraphData data;

      // set thread-appropriate graph data
      GraphDataProvider::exVertglbtab = &data.exVertglbtabRaw[0];
      int vertloccnt = exVertglbtab[mpiMe+1] - exVertglbtab[mpiMe];

      int vertIdx = 0;
      int edgeIdx = 0;
      for (int i=0; i<mpiMe; ++i) {
        vertIdx += exVertglbtab[i+1] - exVertglbtab[i] + 1;
        edgeIdx += data.exVertloctabRaw[vertIdx-1];
      }
      GraphDataProvider::exVertloctab = &data.exVertloctabRaw[vertIdx];
      GraphDataProvider::exVertwgtcmp = &data.exVertwgtcmpRaw[vertIdx];
      GraphDataProvider::exVertwgtvis = &data.exVertwgtvisRaw[vertIdx];

      GraphDataProvider::exEdgeloctab = &data.exEdgeloctabRaw[edgeIdx];
      GraphDataProvider::exEdgewgtcmp = &data.exEdgewgtcmpRaw[edgeIdx];
      GraphDataProvider::exEdgewgtvis = &data.exEdgewgtvisRaw[edgeIdx];

      if (VERBOSE) {
        COUTMPIME << "[ vl v# e# ]" << vertloccnt << " " << vertIdx << " " << edgeIdx << endl;
        COUTMPIME << "[ vloctab  ]";
        for (int i=0; i<=vertloccnt; ++i) { cout << " " << exVertloctab[i]; }
        cout << endl;
        COUTMPIME << "[ vwgt cmp ]";
        for (int i=0; i<=vertloccnt; ++i) { cout << " " << exVertwgtcmp[i]; }
        cout << endl;
        COUTMPIME << "[ vwgt vis ]";
        for (int i=0; i<=vertloccnt; ++i) { cout << " " << exVertwgtvis[i]; }
        cout << endl;
        COUTMPIME << "[ eloctab  ]";
        for (int i=0; i<exVertloctab[vertloccnt]; ++i) { cout << " " << exEdgeloctab[i]; }
        cout << endl;
        COUTMPIME << "[ ewgt cmp ]";
        for (int i=0; i<exVertloctab[vertloccnt]; ++i) { cout << " " << exEdgewgtcmp[i]; }
        cout << endl;
        COUTMPIME << "[ ewgt vis ]";
        for (int i=0; i<exVertloctab[vertloccnt]; ++i) { cout << " " << exEdgewgtvis[i]; }
        cout << endl;
      }
    }

};

// init static variables
int GraphDataProvider::mpiN = 0;
int GraphDataProvider::mpiMe = 0;

int* GraphDataProvider::exVertglbtab = NULL;
int* GraphDataProvider::exVertloctab = NULL;
int* GraphDataProvider::exEdgeloctab = NULL;
int* GraphDataProvider::exVertwgtcmp = NULL;
int* GraphDataProvider::exEdgewgtcmp = NULL;
int* GraphDataProvider::exVertwgtvis = NULL;
int* GraphDataProvider::exEdgewgtvis = NULL;

class MainSystemTest : public testing::Test, protected GraphDataProvider {
  protected:
    PPStee* _ppstee;
    static PPSteeGraph* _pgraph;
    static PPSteeWeights* _pwgtCmp;
    static PPSteeWeights* _pwgtVis;
    static PPSteePart* _ppart;

    static void SetUpTestCase() {

      initGraphData();

      if (VERBOSE) { COUTMPIME << "Starting..." << endl; }

      // get graph (as ParMETIS type)
      _pgraph = new PPSteeGraphParmetis(MPI_COMM_WORLD, exVertglbtab, exVertloctab, exEdgeloctab, PPSTEE_DATAACCESS_COPY);

      if (VERBOSE) { COUTMPIME << "Graph created." << endl; }

      // construct weights for computation
      _pwgtCmp = new PPSteeWeights(_pgraph);
      _pwgtCmp->setWeightsData(exVertwgtcmp, exEdgewgtcmp);

      if (VERBOSE) { COUTMPIME << "CmpWeights created." << endl; }

      // construct weights for visualisation
      _pwgtVis = new PPSteeWeights(_pgraph);
      _pwgtVis->setWeightsData(exVertwgtvis, exEdgewgtvis);

      if (VERBOSE) { COUTMPIME << "VisWeights created." << endl; }

      // construct partitioning
      _ppart = new PPSteePart(_pgraph);

      if (VERBOSE) { COUTMPIME << "Part created." << endl; }

    };

    static void TearDownTestCase() {
      delete _pgraph;
      delete _pwgtCmp;
      delete _pwgtVis;
      delete _ppart;
    }

    virtual void SetUp() {

      // get interface
      _ppstee = new PPStee();

      if (VERBOSE) { COUTMPIME << "PPSteeObj created." << endl; }

    };

    virtual void TearDown() {
      delete _ppstee;
    }

  void printGraphData(ostream& os, PPSteeGraph* pgraph) const {
    OSMPIME;
    os << "vertloccnt: " << pgraph->getVertloccnt();
    os << " edgeloccnt: " << pgraph->getEdgeloccnt();
    os << endl;
    OSMPIME << "vertloctab:";
    for (int i=0; i<=pgraph->getVertloccnt(); ++i) { os << " " << pgraph->getVertloctab()[i]; }
    os << endl;
    OSMPIME << "edgeloctab:";
    for (int i=0; i<pgraph->getEdgeloccnt(); ++i) { os << " " << pgraph->getEdgeloctab()[i]; }
    os << endl;
  }

};

// init static variables
PPSteeGraph* MainSystemTest::_pgraph = NULL;
PPSteeWeights* MainSystemTest::_pwgtCmp = NULL;
PPSteeWeights* MainSystemTest::_pwgtVis = NULL;
PPSteePart* MainSystemTest::_ppart = NULL;


// Tests

TEST_F(MainSystemTest, SystemTestPartitionerParmetisGraphParmetis) {
  // gather threads
  MPI_Barrier(MPI_COMM_WORLD);

  // submit graph
  _ppstee->submitGraph(*_pgraph);

  if (VERBOSE) { COUTMPIME << "Graph submitted." << endl; }

  // submit weights
  _ppstee->submitNewStage(*_pwgtCmp, PPSTEE_STAGE_COMPUTATION);
  _ppstee->submitNewStage(*_pwgtVis, PPSTEE_STAGE_VISUALISATION);

  if (VERBOSE) { COUTMPIME << "Weights submitted." << endl; }

  // get partitioning
  _ppstee->getPartitioning(&_ppart);

  // check partitioning
  ASSERT_TRUE(_ppart);
  ASSERT_TRUE(_ppart->getPartData());
  for (int i=0; i<_ppart->getVertloccnt(); ++i) {
    ASSERT_LE(0, _ppart->getPartData()[i]);
    ASSERT_LT(_ppart->getPartData()[i], mpiN);
  }
}

TEST_F(MainSystemTest, SystemTestPartitionerPtscotchGraphPtscotchConvertedFromParmetis) {
  // gather threads
  MPI_Barrier(MPI_COMM_WORLD);

  // submit graph
  _ppstee->submitGraph(*_pgraph);

  if (VERBOSE) { COUTMPIME << "Graph submitted." << endl; }

  // submit weights
  _ppstee->submitNewStage(*_pwgtCmp, PPSTEE_STAGE_COMPUTATION);
  _ppstee->submitNewStage(*_pwgtVis, PPSTEE_STAGE_VISUALISATION);

  if (VERBOSE) { COUTMPIME << "Weights submitted." << endl; }

  // convert graph to PTScotch
  PPSteeGraph newGraph(_pgraph, PPSTEE_GRAPH_PTSCOTCH);

  printGraphData(cout, &newGraph);
  cout << endl;

  // resubmit new graph
  _ppstee->submitGraph(newGraph);

  // calculate new partitioning
  PPSteePart* newPart = NULL;
  _ppstee->getPartitioning(&newPart);

  // check partitioning
  ASSERT_TRUE(newPart);
  ASSERT_TRUE(newPart->getPartData());
  int mpiMe = newGraph.getMpiMe();
  COUTMPIME << "[ part sco ]";
  for (int i=0; i<newPart->getVertloccnt(); ++i) {
    cout << " " << newPart->getPartData()[i];
    ASSERT_LE(0, newPart->getPartData()[i]);
    ASSERT_LT(newPart->getPartData()[i], mpiN);
  }
  cout << endl;
}

TEST_F(MainSystemTest, SystemTestPartitionerZoltanGraphZoltanConvertedFromParmetis) {
  // gather threads
  MPI_Barrier(MPI_COMM_WORLD);

  // submit graph
  _ppstee->submitGraph(*_pgraph);

  if (VERBOSE) { COUTMPIME << "Graph submitted." << endl; }

  // submit weights
  _ppstee->submitNewStage(*_pwgtCmp, PPSTEE_STAGE_COMPUTATION);
  _ppstee->submitNewStage(*_pwgtVis, PPSTEE_STAGE_VISUALISATION);

  if (VERBOSE) { COUTMPIME << "Weights submitted." << endl; }

  // convert graph to Zoltan
  PPSteeGraph newGraphZoltan(_pgraph, PPSTEE_GRAPH_ZOLTAN);

  //  printGraphData(cout, &newGraphZoltan);
  cout << endl;

  // resubmit new graph
  _ppstee->submitGraph(newGraphZoltan);

  // calculate new partitioning
  PPSteePart* newPartZoltan = NULL;
  _ppstee->getPartitioning(&newPartZoltan);

  // check partitioning
  ASSERT_TRUE(newPartZoltan);
  ASSERT_TRUE(newPartZoltan->getPartData());
  int mpiMe = newGraphZoltan.getMpiMe();
  COUTMPIME << "[ part zol ]";
  for (int i=0; i<newPartZoltan->getVertloccnt(); ++i) {
    cout << " " << newPartZoltan->getPartData()[i];
    ASSERT_LE(0, newPartZoltan->getPartData()[i]);
    ASSERT_LT(newPartZoltan->getPartData()[i], mpiN);
  }
  cout << endl;
}
