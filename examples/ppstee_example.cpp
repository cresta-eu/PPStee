/*!
 * \file
 *
 * \brief PPStee example: main usage
 * \details This example demonstrates basic usage of PPStee.
 * The example graph was generated with graphGen (see tools/graphGen) and is similar to the graph used in the ParMETIS manual (cf. http://glaros.dtc.umn.edu/gkhome/fetch/sw/parmetis/manual.pdf);
 * basically it is a 7x9 version of the 3x5 rectangle-shaped regular mesh with row-major numbering of the vertices used in the ParMETIS manual.
 *
 * \note This example MUST be run on 3 MPI threads.
 *
 * \author Gregor Matura
 *
 */

#include "ppstee.hpp"

#include <iostream>
#include <unistd.h>


#define VERBOSE 1
#ifndef GDBDEBUG
  #define GDBDEBUG 0
#endif /* GDBDEBUG */

#define COUTMPIME cout << "[" << mpiMe << "] "

using std::cout;
using std::endl;
using std::ostream;


// include graph data
#include "test/res/graphData-3.cpp"


int main(int argc, char** argv) {

  // MPI init
  MPI_Init(&argc, &argv);
  int mpiMe, mpiN;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpiMe);
  MPI_Comm_size(MPI_COMM_WORLD, &mpiN);
  if (mpiN!=3) {
    cout << "Start with 3 MPI threads! Abort." << endl;
    return -1;
  }

#if GDBDEBUG != 0
  // prepare on-start breakpoint
  int i=0;
  cout << "[" << mpiMe << ":pidinfo] " << getpid() << endl;
  cout.flush();
  while (i == 0) {
    sleep(5);
  }
#endif /* GDBDEBUG != 0 */

  // init example graph data
  int* exVertglbtab;
  int* exVertloctab;
  int* exEdgeloctab;
  int* exVertwgtcmp;
  int* exEdgewgtcmp;
  int* exVertwgtvis;
  int* exEdgewgtvis;
  
  // get example graph data
  exGraphData data;

  // set thread-appropriate graph data
  exVertglbtab = &data.exVertglbtabRaw[0];
  int vertloccnt = exVertglbtab[mpiMe+1] - exVertglbtab[mpiMe];

  int vertIdx = 0;
  int edgeIdx = 0;
  for (int i=0; i<mpiMe; ++i) {
    vertIdx += exVertglbtab[i+1] - exVertglbtab[i] + 1;
    edgeIdx += data.exVertloctabRaw[vertIdx-1];
  }
  exVertloctab = &data.exVertloctabRaw[vertIdx];
  exVertwgtcmp = &data.exVertwgtcmpRaw[vertIdx];
  exVertwgtvis = &data.exVertwgtvisRaw[vertIdx];

  exEdgeloctab = &data.exEdgeloctabRaw[edgeIdx];
  exEdgewgtcmp = &data.exEdgewgtcmpRaw[edgeIdx];
  exEdgewgtvis = &data.exEdgewgtvisRaw[edgeIdx];

  if (VERBOSE) {
    COUTMPIME << vertloccnt << " " << vertIdx << " " << edgeIdx << endl;
    COUTMPIME;
    for (int i=0; i<=vertloccnt; ++i) { cout << " " << exVertloctab[i]; }
    cout << endl;
    COUTMPIME;
    for (int i=0; i<=vertloccnt; ++i) { cout << " " << exVertwgtcmp[i]; }
    cout << endl;
    COUTMPIME;
    for (int i=0; i<=vertloccnt; ++i) { cout << " " << exVertwgtvis[i]; }
    cout << endl;
    COUTMPIME;
    for (int i=0; i<exVertloctab[vertloccnt]; ++i) { cout << " " << exEdgeloctab[i]; }
    cout << endl;
    COUTMPIME;
    for (int i=0; i<exVertloctab[vertloccnt]; ++i) { cout << " " << exEdgewgtcmp[i]; }
    cout << endl;
    COUTMPIME;
    for (int i=0; i<exVertloctab[vertloccnt]; ++i) { cout << " " << exEdgewgtvis[i]; }
    cout << endl;
  }

  // get graph (as ParMETIS type)
  if (VERBOSE) { COUTMPIME << "Starting..." << endl; }

  // get graph (as ParMETIS type)
  PPSteeGraph* _pgraph = new PPSteeGraphParmetis(MPI_COMM_WORLD, exVertglbtab, exVertloctab, exEdgeloctab, PPSTEE_DATAACCESS_COPY);

  if (VERBOSE) { COUTMPIME << "Graph created." << endl; }

  // construct and set weights for computation
  PPSteeWeights* _pwgtCmp = new PPSteeWeights(_pgraph);
  _pwgtCmp->setWeightsData(exVertwgtcmp, exEdgewgtcmp);

  if (VERBOSE) { COUTMPIME << "CmpWeights created." << endl; }

  // construct and set weights for visualisation
  PPSteeWeights* _pwgtVis = new PPSteeWeights(_pgraph);
  _pwgtVis->setWeightsData(exVertwgtvis, exEdgewgtvis);

  if (VERBOSE) { COUTMPIME << "VisWeights created." << endl; }

  // get interface
  PPStee* _ppstee = new PPStee();

  if (VERBOSE) { COUTMPIME << "PPSteeObj created." << endl; }

  // submit graph
  _ppstee->submitGraph(*_pgraph);

  if (VERBOSE) { COUTMPIME << "Graph submitted." << endl; }

  // submit weights
  _ppstee->submitNewStage(*_pwgtCmp, PPSTEE_STAGE_COMPUTATION);
  _ppstee->submitNewStage(*_pwgtVis, PPSTEE_STAGE_VISUALISATION);

  if (VERBOSE) { COUTMPIME << "Weights submitted." << endl; }

  // calculate partitioning
  PPSteePart* _ppart = new PPSteePart(_pgraph);
  _ppstee->getPartitioning(&_ppart);

  if (VERBOSE) { COUTMPIME << "Partitioning created with ParMETIS." << endl; }

  // output
  COUTMPIME << "vertex redist. by ParMETIS: ";
  if (_ppart && _ppart->getPartData()) for (int i=0; i<_ppart->getVertloccnt(); ++i) cout << (i==0?"":" | ") << i << " -> " << _ppart->getPartData()[i];
  cout << endl;

  // convert graph to PTScotch
  PPSteeGraph newGraph(_pgraph, PPSTEE_GRAPH_PTSCOTCH);

  if (VERBOSE) { COUTMPIME << "Graph converted to PTScotch type." << endl; }

  // resubmit new graph
  _ppstee->submitGraph(newGraph);

  if (VERBOSE) { COUTMPIME << "New graph resubmitted." << endl; }

  // calculate new partitioning
  PPSteePart* newPart = NULL;
  _ppstee->getPartitioning(&newPart);

  if (VERBOSE) { COUTMPIME << "Partitioning created with PTScotch." << endl; }

  // output
  COUTMPIME << "vertex redist. by PTScotch: ";
  if (newPart!=NULL && newPart->getPartData()) for (int i=0; i<newPart->getVertloccnt(); ++i) cout << (i==0?"":" | ") << i << " -> " << newPart->getPartData()[i];
  cout << endl;

  // convert graph to Zoltan
  PPSteeGraph newGraphZoltan(_pgraph, PPSTEE_GRAPH_ZOLTAN);

  if (VERBOSE) { COUTMPIME << "Graph converted to Zoltan type." << endl; }

  // resubmit new graph
  _ppstee->submitGraph(newGraphZoltan);

  if (VERBOSE) { COUTMPIME << "New graph resubmitted." << endl; }

  // calculate new partitioning
  PPSteePart* newPartZoltan = NULL;
  _ppstee->getPartitioning(&newPartZoltan);

  if (VERBOSE) { COUTMPIME << "Partitioning created with Zoltan." << endl; }

  // output
  COUTMPIME << "vertex redist. by Zoltan: ";
  if (newPartZoltan!=NULL && newPartZoltan->getPartData()) for (int i=0; i<newPartZoltan->getVertloccnt(); ++i) cout << (i==0?"":" | ") << i << " -> " << newPartZoltan->getPartData()[i];
  cout << endl;

  if (VERBOSE) { COUTMPIME << "Done." << endl; }

  // MPI finalize
  MPI_Finalize();

  // exit
  return 0;
}
