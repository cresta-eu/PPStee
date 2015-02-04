/*! \file ppstee_objects.cpp
 * \brief Defining objects used by PPStee.
 *
 * \date 30.01.2013
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


#include "ppstee_objects.hpp"


PPSteeGraph::PPSteeGraph(MPI_Comm mpiComm, PPSTEE_DATAACCESS access/* = PPSTEE_DATAACCESS_COPY */)
  : type(PPSTEE_GRAPH_UNDETERMINED),
    dataAccess(access),
    onlyLocalData(true),
    mpiComm(mpiComm),
    vertglbtab(NULL),
    vertloccnt(0),
    edgeloccnt(0),
    vertloctab(NULL),
    edgeloctab(NULL) {

  // init MPI values
  MPI_Comm_rank(mpiComm, &mpiMe);
  MPI_Comm_size(mpiComm, &mpiN);
}

PPSteeGraph::PPSteeGraph(PPSteeGraph* existingGraph, PPSTEE_GRAPH targetType)
: type(targetType),
  dataAccess(existingGraph->getDataAccess()),
  onlyLocalData(existingGraph->getOnlyLocalData()),
  mpiComm(existingGraph->getMpiComm()),
  mpiMe(existingGraph->getMpiMe()),
  mpiN(existingGraph->getMpiN()),
  vertglbtab(NULL),
  vertloccnt(existingGraph->getVertloccnt()),
  edgeloccnt(existingGraph->getEdgeloccnt()),
  vertloctab(NULL),
  edgeloctab(NULL) {

  // if existing graph has Zoltan type \em and uses query functions (#dataAccess) set new #dataAccess to #PPSTEE_DATAACCESS_COPY;
  // due to the upcoming use of getVertloctab() a data copy operation will be forced.
  if (existingGraph->getType() == PPSTEE_GRAPH_ZOLTAN &&
      existingGraph->getDataAccess() == PPSTEE_DATAACCESS_QRYFN) {
    dataAccess = PPSTEE_DATAACCESS_COPY;
  }

  // copy/reference array data
  if (onlyLocalData) {
    PPSTEE_ERROR retError = setLocalGraphData(existingGraph->getVertloccnt(), existingGraph->getVertloctab(), existingGraph->getEdgeloctab());
  } else {
    PPSTEE_ERROR retError = setGraphData(existingGraph->getVertglbtab(), existingGraph->getVertloctab(), existingGraph->getEdgeloctab());
  }
}

PPSteeGraph::~PPSteeGraph() {
  freeMemoryIfHoldingDataCopies();
}

PPSTEE_ERROR PPSteeGraph::refineMeshAutomatically() {return PPSTEE_ERROR_UNKNOWNERROR; }  // TODO: implement


PPSTEE_ERROR PPSteeGraph::setGraphData(int* vertglbtab, int* vertloctab, int* edgeloctab) {
  // set reference to or copy data (first for global vertex table, then for local vertex and edge tables)
  if (dataAccess==PPSTEE_DATAACCESS_VIEW) {
    this->vertglbtab = vertglbtab;
  } else {
    int vertglbsize = 4*(mpiN+1);
    this->vertglbtab = (int*)malloc(vertglbsize);
    memcpy(this->vertglbtab, vertglbtab, vertglbsize);
  }
  setLocalGraphData(vertglbtab[mpiMe+1] - vertglbtab[mpiMe], vertloctab, edgeloctab);

  // set local modifier (has to be at the end of this function; setLocalGraphData() also modifies this attribute.)
  onlyLocalData = false;

  return PPSTEE_ERROR_SUCCESS;  // TODO: error handling
}

PPSTEE_ERROR PPSteeGraph::setLocalGraphData(int vertloccnt, int* vertloctab, int* edgeloctab) {
  // set local modifier
  onlyLocalData = true;

  // set reference to or copy data (only for loca vertex and edge tables)
  this->vertloccnt = vertloccnt;
  this->edgeloccnt = vertloctab[vertloccnt];
  if (dataAccess==PPSTEE_DATAACCESS_VIEW) {
    this->vertloctab = vertloctab;
    this->edgeloctab = edgeloctab;
  } else {
    int vertlocsize = 4*(vertloccnt+1);
    this->vertloctab = (int*)malloc(vertlocsize);
    memcpy(this->vertloctab, vertloctab, vertlocsize);
    int edgelocsize = 4*edgeloccnt;
    this->edgeloctab = (int*)malloc(edgelocsize);
    memcpy(this->edgeloctab, edgeloctab, edgelocsize);
  }

  return PPSTEE_ERROR_SUCCESS;  // TODO: error handling
}

PPSTEE_ERROR PPSteeGraph::buildGlobalDataFromLocalData() {
  // check if necessary or impossible
  if (!onlyLocalData || vertloccnt<=0) {
    return PPSTEE_ERROR_SUCCESS;
  }

  // malloc mem
  if (vertglbtab != NULL) {
    free (vertglbtab);
  }
  vertglbtab = (int*) malloc((mpiN+1)*4);

  // gather all #vertloccnt
  MPI_Allgather(&vertloccnt, 1, MPI_INT, &(vertglbtab[1]), 1, MPI_INT, mpiComm);

  // calculate correct #vertglbtab
  vertglbtab[0] = 0;
  for (int i=1; i<mpiN; ++i) {
    vertglbtab[i+1] += vertglbtab[i];
  }

  return PPSTEE_ERROR_SUCCESS;
}

void PPSteeGraph::setType(PPSTEE_GRAPH type) {
  this->type = type;
}

int* PPSteeGraph::getVertglbtab() {
  if (onlyLocalData) {
    buildGlobalDataFromLocalData();
  }
  return vertglbtab;
}

void PPSteeGraph::freeMemoryIfHoldingDataCopies() {
  // free memory if holding data copies
  if (dataAccess==PPSTEE_DATAACCESS_COPY) {
    if (vertglbtab!=NULL) { free(vertglbtab); }
    if (vertloctab!=NULL) { free(vertloctab); }
    if (edgeloctab!=NULL) { free(edgeloctab); }
  }
}

PPSteeGraphParmetis::PPSteeGraphParmetis(MPI_Comm mpiComm, int* vtxdist, int* xadj, int* adjncy, PPSTEE_DATAACCESS access /* = PPSTEE_DATAACCESS_VIEW */)
  : PPSteeGraph(mpiComm, access) {

  create(vtxdist, xadj, adjncy);
  setType(PPSTEE_GRAPH_PARMETIS);
}

PPSteeGraphParmetis::PPSteeGraphParmetis(PPSteeGraph* existingGraph)
  : PPSteeGraph(existingGraph, PPSTEE_GRAPH_PARMETIS) {};

PPSTEE_ERROR PPSteeGraphParmetis::create(int* vtxdist, int* xadj, int* adjncy) {
  // arrange data so that it can be submitted
  // for ParMETIS: (almost) nothing to do
  return setGraphData(vtxdist, xadj, adjncy);
}

PPSteeGraphPtscotch::PPSteeGraphPtscotch(MPI_Comm mpiComm,
                    int  baseval,
                    int  vertlocnbr,
                    int  vertlocmax,
                    int* vertloctab,
                    int* vendloctab,
                    int* veloloctab,
                    int* vlblloctab,
                    int  edgelocnbr,
                    int  edgelocsiz,
                    int* edgeloctab,
                    int* edgegsttab,
                    int* edloloctab,
                    PPSTEE_DATAACCESS access/* = PPSTEE_DATAACCESS_COPY */
                    )
  : PPSteeGraph(mpiComm, access) {
  create(vertlocnbr, vertloctab, edgeloctab);
  setType(PPSTEE_GRAPH_PTSCOTCH);
}

PPSteeGraphPtscotch::PPSteeGraphPtscotch(MPI_Comm mpiComm,
                    int vertlocnbr,
                    int* vertloctab,
                    int* edgeloctab,
                    PPSTEE_DATAACCESS access/* = PPSTEE_DATAACCESS_COPY */
                    )
  : PPSteeGraph(mpiComm, access) {
  create(vertlocnbr, vertloctab, edgeloctab);
  setType(PPSTEE_GRAPH_PTSCOTCH);
}

PPSteeGraphPtscotch::PPSteeGraphPtscotch(PPSteeGraph* existingGraph)
: PPSteeGraph(existingGraph, PPSTEE_GRAPH_PTSCOTCH) {};

PPSTEE_ERROR PPSteeGraphPtscotch::create(int vertlocnbr, int* vertloctab, int* edgeloctab) {
  // arrange data so that it can be submitted
  // for PTScotch: (almost) nothing to do
  return setLocalGraphData(vertlocnbr, vertloctab, edgeloctab);
}

PPSteeGraphZoltan::PPSteeGraphZoltan(MPI_Comm mpiComm,
      ZOLTAN_NUM_OBJ_FN * zNumObjFn,
      ZOLTAN_OBJ_LIST_FN* zObjListFn,
      ZOLTAN_NUM_EDGES_MULTI_FN* zNumEdgesMultiFn,
      ZOLTAN_EDGE_LIST_MULTI_FN* zEdgeListMultiFn,
      void* dNumObj/* = NULL*/,
      void* dObjList/* = NULL */,
      void* dNumEdgesMulti/* = NULL */,
      void* dEdgeListMulti/* = NULL */,
      PPSTEE_DATAACCESS access/* = PPSTEE_DATAACCESS_COPY */
      )
  : PPSteeGraph(mpiComm, access) {
  create(zNumObjFn, zObjListFn, zNumEdgesMultiFn, zEdgeListMultiFn);
  setType(PPSTEE_GRAPH_ZOLTAN);
}


PPSteeGraphZoltan::PPSteeGraphZoltan(PPSteeGraph* existingGraph)
: PPSteeGraph(existingGraph, PPSTEE_GRAPH_ZOLTAN) {};

/**
 *  Arrange data so that it can be submitted;
 *  for Zoltan:
 *  if #dataAccess == #PPSTEE_DATAACCESS_COPY: create vertglbtab, vertloctab and edgeloctab using provided Zoltan query functions.
 *  if #dataAccess == #PPSTEE_DATAACCESS_QRYFN: save references to Zoltan query functions for later use.
 *  if #dataAccess == #PPSTEE_DATAACCESS_VIEW: same as #PPSTEE_DATAACCESS_QRYFN.
 *
 *  Additionally set #vertloccnt and #edgeloccnt.
 *
 */
PPSTEE_ERROR PPSteeGraphZoltan::create(
      ZOLTAN_NUM_OBJ_FN * zNumObjFn,
      ZOLTAN_OBJ_LIST_FN* zObjListFn,
      ZOLTAN_NUM_EDGES_MULTI_FN* zNumEdgesMultiFn,
      ZOLTAN_EDGE_LIST_MULTI_FN* zEdgeListMultiFn,
      void* dNumObj/* = NULL*/,
      void* dObjList/* = NULL */,
      void* dNumEdgesMulti/* = NULL */,
      void* dEdgeListMulti/* = NULL */){

  this->zNumObjFn = zNumObjFn;
  this->zObjListFn = zObjListFn;
  this->zNumEdgesMultiFn = zNumEdgesMultiFn;
  this->zEdgeListMultiFn = zEdgeListMultiFn;

  this->dNumObj = dNumObj;
  this->dObjList = dObjList;
  this->dNumEdgesMulti = dNumEdgesMulti;
  this->dEdgeListMulti = dEdgeListMulti;


  if (dataAccess == PPSTEE_DATAACCESS_VIEW) {
    dataAccess = PPSTEE_DATAACCESS_QRYFN;
  } else if (dataAccess == PPSTEE_DATAACCESS_COPY) {
      return createDataCopiesFromQueryFunctions();
  }

  return PPSTEE_ERROR_SUCCESS;
}

PPSTEE_ERROR PPSteeGraphZoltan::createDataCopiesFromQueryFunctions() {
  return queryZoltanFunctionsForLocalGraphData();
}

PPSTEE_ERROR PPSteeGraphZoltan::queryZoltanFunctionsForLocalGraphData() {
  int ierr;

  // get #vertloccnt
  vertloccnt = zNumObjFn(dNumObj, &ierr);
  if (ierr!=ZOLTAN_OK) { return PPSTEE_ERROR_ZOLTAN; }

  // prepare for querying
  ZOLTAN_ID_PTR global_ids;
  ZOLTAN_ID_PTR local_ids;
  global_ids = (ZOLTAN_ID_PTR) malloc(vertloccnt*4);
  local_ids = (ZOLTAN_ID_PTR) malloc(vertloccnt*4);

  // use query function #zObjListFn to fill global_ids and local_ids
  zObjListFn(dObjList, PPSTEE_NUM_GID_ENTRIES, PPSTEE_NUM_LID_ENTRIES, global_ids, local_ids, 0, NULL, &ierr);

  // query #zNumEdgesMultiFn for num_edges, use this for #vertloctab and #edgeloccnt
  int* num_edges = (int*) malloc(vertloccnt*4);
  zNumEdgesMultiFn(dNumEdgesMulti, PPSTEE_NUM_GID_ENTRIES, PPSTEE_NUM_LID_ENTRIES, vertloccnt, global_ids, local_ids, num_edges, &ierr);

  // assemble #vertloctab and #edgeloccnt
  vertloctab = (int*) malloc((vertloccnt+1)*4);
  vertloctab[0] = 0;
  edgeloccnt = 0;
  for (int i=0; i<vertloccnt; ++i) {
    vertloctab[i+1] += vertloctab[i] + num_edges[i];
    edgeloccnt += num_edges[i];
  }

  // query #zEdgeListMultiFn for #edgeloctab
  edgeloctab = (int*) malloc(edgeloccnt*4);
  ZOLTAN_ID_PTR nbor_global_id = (ZOLTAN_ID_PTR) edgeloctab;
  int* nbor_procs = (int*) malloc(edgeloccnt*4);
  zEdgeListMultiFn(dEdgeListMulti, PPSTEE_NUM_GID_ENTRIES, PPSTEE_NUM_LID_ENTRIES, vertloccnt, global_ids, local_ids, num_edges, nbor_global_id, nbor_procs, 0, NULL, &ierr);

  // free mem
  free(global_ids);
  free(local_ids);
  free(num_edges);
  free(nbor_procs);

  // we now have graph data, but only locally
  dataAccess = PPSTEE_DATAACCESS_COPY;
  onlyLocalData = true;

  return (ierr==ZOLTAN_OK ? PPSTEE_ERROR_SUCCESS : PPSTEE_ERROR_ZOLTAN );
}

int* PPSteeGraphZoltan::getVertglbtab() {
  if (dataAccess==PPSTEE_DATAACCESS_QRYFN) {
      createDataCopiesFromQueryFunctions();
      buildGlobalDataFromLocalData();
  }
  return vertglbtab;
}

int* PPSteeGraphZoltan::getVertloctab() {
  if (dataAccess==PPSTEE_DATAACCESS_QRYFN) {
      createDataCopiesFromQueryFunctions();
  }
  return vertloctab;
}

int* PPSteeGraphZoltan::getEdgeloctab() {
  if (dataAccess==PPSTEE_DATAACCESS_QRYFN) {
      createDataCopiesFromQueryFunctions();
  }
  return edgeloctab;
}

PPSteeWeights::PPSteeWeights(PPSteeGraph* graph)
  : type(PPSTEE_WEIGHTS_ALL),
    dataAccess(graph->getDataAccess()),
    vertloccnt(graph->getVertloccnt()),
    edgeloccnt(graph->getEdgeloccnt()),
    vertexWeights(NULL),
    edgeWeights(NULL),
    vertexWeightsEquipartitioned(true),
    edgeWeightsEquipartitioned(true) {
}

PPSteeWeights::~PPSteeWeights() {
  freeMemoryIfHoldingDataCopies();
}

PPSTEE_ERROR PPSteeWeights::setWeightsData(int* vertexWeights, int* edgeWeights) {
  PPSTEE_ERROR retErrVertex, retErrEdge;

  if (type == PPSTEE_WEIGHTS_ALL || type == PPSTEE_WEIGHTS_ONLYVERTEX) {
    retErrVertex = setOnlyVertexWeightsData(vertexWeights);
  }

  if (type == PPSTEE_WEIGHTS_ALL || type == PPSTEE_WEIGHTS_ONLYEDGE) {
    retErrEdge = setOnlyEdgeWeightsData(edgeWeights);
  }

  return (retErrVertex==PPSTEE_ERROR_SUCCESS && retErrEdge==PPSTEE_ERROR_SUCCESS ? PPSTEE_ERROR_SUCCESS : PPSTEE_ERROR_UNKNOWNERROR);  // TODO: error handling
}

PPSTEE_ERROR PPSteeWeights::setOnlyVertexWeightsData(int* vertexWeights) {
  if (dataAccess==PPSTEE_DATAACCESS_VIEW) {
    this->vertexWeights = vertexWeights;
  } else {
    if (this->vertexWeights!=NULL) { free(this->vertexWeights); }
    this->vertexWeights = (int*) malloc(4*vertloccnt);
    memcpy(this->vertexWeights, vertexWeights, 4*vertloccnt);
  }

  vertexWeightsEquipartitioned = false;

  return PPSTEE_ERROR_SUCCESS;  // TODO: error handling
}

PPSTEE_ERROR PPSteeWeights::setOnlyEdgeWeightsData(int* edgeWeights) {
  if (dataAccess==PPSTEE_DATAACCESS_VIEW) {
    this->edgeWeights = edgeWeights;
  } else {
    if (this->edgeWeights!=NULL) { free(this->edgeWeights); }
    this->edgeWeights = (int*) malloc(4*edgeloccnt);
    memcpy(this->edgeWeights, edgeWeights, 4*edgeloccnt);
  }

  edgeWeightsEquipartitioned = false;

  return PPSTEE_ERROR_SUCCESS;  // TODO: error handling
}

void PPSteeWeights::freeMemoryIfHoldingDataCopies() {
  // free memory if holding data copies
  if (dataAccess==PPSTEE_DATAACCESS_COPY) {
    if (vertexWeights!=NULL) { free(vertexWeights); }
    if (edgeWeights!=NULL) { free(edgeWeights); }
  }
}

PPSteePart::PPSteePart(PPSteeGraph* graph, int* existingPartitioning)
  : vertloccnt(graph->getVertloccnt()), partitioning(NULL) {
  PPSTEE_ERROR ret = setPartitioningData(existingPartitioning);
}

PPSteePart::PPSteePart(PPSteeGraph* graph)
  : vertloccnt(graph->getVertloccnt()), partitioning(NULL) {
}

//! Destructor
PPSteePart::~PPSteePart() {
  freeMemoryIfHoldingDataCopies();
}


PPSTEE_ERROR PPSteePart::setPartitioningData(int* partData) {
  if (partitioning!=NULL) { free(partitioning); }
  partitioning = (int*) malloc(4*vertloccnt);
  memcpy(partitioning, partData, 4*vertloccnt);

  return PPSTEE_ERROR_SUCCESS;  // TODO: error handling
}

void PPSteePart::freeMemoryIfHoldingDataCopies() {
  // free memory if holding data copies
  if (partitioning!=NULL) { free(partitioning); }
}
