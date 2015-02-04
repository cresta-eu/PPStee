/*! \file ppstee.cpp
 * \brief Main source file.
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

#include "ppstee.hpp"

PPStee::PPStee() :
    graph(NULL), useOldPartitioning(false), oldPart(NULL) {
}

PPSTEE_ERROR PPStee::submitGraph(PPSteeGraph& graph) {
  if (this->graph == NULL) {
    // graph is empty, set to new
    this->graph = &graph;
  } else {
    // a graph is present: check graph-weight compatibility
    PPSteeGraph* oldGraph = this->graph;
    this->graph = &graph;
    if (!isCompatibleGraphAndWeights()) {
      // graph and weights not compatible: restore old graph pointer and return error
      this->graph = oldGraph;
      return PPSTEE_ERROR_MATCHINGERROR;
    }
  }

  return PPSTEE_ERROR_SUCCESS; // TODO: error handling
}

PPSTEE_ERROR PPStee::submitNewStageByWeights(PPSteeWeights& weights, PPSTEE_STAGE type, int* stageNumber /* = NULL*/) {
  // add weights to list
  stages.push_back(pair<PPSTEE_STAGE, PPSteeWeights*>(type, &weights));

  // return stageNumber if requested
  if (stageNumber != NULL) {
    *stageNumber = stages.size() - 1;
  }

  return PPSTEE_ERROR_SUCCESS; // TODO: error handling
}

PPSTEE_ERROR PPStee::changeStageByWeights(PPSteeWeights& weights,
    int stageNumber) {
  // change weights
  stages[stageNumber].second = &weights;

  return PPSTEE_ERROR_SUCCESS; // TODO: error handling
}

PPSTEE_ERROR PPStee::removeStage(int stageNumber) {
  // remove stage: set PPSTEE_STAGE to removed, free mem, reset pointer
  stages[stageNumber].first = PPSTEE_STAGE_REMOVED;
  stages[stageNumber].second->freeMemoryIfHoldingDataCopies();
  stages[stageNumber].second = NULL;

  return PPSTEE_ERROR_SUCCESS; // TODO: error handling
}

PPSTEE_ERROR PPStee::submitOldPartitioning(PPSteePart& part) {
  this->oldPart = &part;

  return PPSTEE_ERROR_SUCCESS; // TODO: error handling
}

int PPStee::getNumberOfActiveStages() const {
  // get all stages and subtract removed
  int nos = getNumberOfAllStages();
  for (vector<pair<PPSTEE_STAGE, PPSteeWeights*> >::const_iterator it = stages.begin(); it != stages.end(); ++it) {
    if (it->first == PPSTEE_STAGE_REMOVED) {
      nos--;
    }
  }
  return nos;
}

PPSTEE_ERROR PPStee::getPartitioning(PPSteePart** part) {
  switch (graph->getType()) {
  case PPSTEE_GRAPH_PARMETIS:
    return getPartitioningParmetis(part);
  case PPSTEE_GRAPH_PTSCOTCH:
    return getPartitioningPtscotch(part);
  case PPSTEE_GRAPH_ZOLTAN:
    return getPartitioningZoltan(part);
  case PPSTEE_GRAPH_UNDETERMINED:
    return PPSTEE_ERROR_GRAPHUNDETERMINED;
  default:
    return PPSTEE_ERROR_UNKNOWNERROR;
  }
}

template<typename T>
PPSTEE_ERROR PPStee::transpose(const int h, const int w, T** array) {
  T* tar = new T[w*h];

  for (int iB = 0; iB < h; iB += PPSTEE_TRANSPOSE_BLOCKSIZE) {
    int iMax = iB + PPSTEE_TRANSPOSE_BLOCKSIZE < h ? iB + PPSTEE_TRANSPOSE_BLOCKSIZE : h;
    for (int jB = 0; jB < w; jB += PPSTEE_TRANSPOSE_BLOCKSIZE) {
      int jMax = jB + PPSTEE_TRANSPOSE_BLOCKSIZE < w ? jB + PPSTEE_TRANSPOSE_BLOCKSIZE : w;
      for (int i = iB; i < iMax; ++i) {
        for (int j = jB; j < jMax; ++j) {
          tar[j*h+i] = (*array)[i*w+j];
        }
      }
    }
  }

  T* tmp = *array;
  *array = tar;
  delete [] tmp;

  return PPSTEE_ERROR_SUCCESS;
}

template<typename T>
PPSTEE_ERROR PPStee::gatherWeightsData(PPSTEE_GRAPHDATUM type, const set<int>& lstStages, const bool accumulate, const bool doAlloc, const int numSpareStages, T** array) {
  int lenLine;
  if (type == PPSTEE_GRAPHDATUM_VERTEX) {
    lenLine = graph->getVertloccnt();
  } else if (type == PPSTEE_GRAPHDATUM_EDGE) {
    lenLine = graph->getEdgeloccnt();
  }
  int lenBlock = lenLine;
  if (not accumulate) { lenBlock *= lstStages.size() + numSpareStages; }

  // alloc
  if (doAlloc) {
    *array = new T[lenBlock];
  }

  // collect explicitly-given weights
  int curIdx = 0;
  for (set<int>::const_iterator it = lstStages.begin(); it != lstStages.end(); ++it) {
    int* wgts = stages[*it].second->getWeights(type);
    if (it == lstStages.begin()) {
      copy(wgts, wgts + lenLine, *array);
    } else {
      if (accumulate) {
        transform(wgts, wgts + lenLine, *array, *array, plus<T>());
      } else {
        copy(wgts, wgts + lenLine, *array + curIdx * lenLine);
      }
    }
    ++curIdx;
  }

  return PPSTEE_ERROR_SUCCESS;
}

PPSTEE_WEIGHTS PPStee::getWeightsFlagAndStageSets(PPSTEE_WEIGHTS targetWgtflag, set<int>& stagesContainingVertexWeights, set<int>& stagesContainingEdgeWeights) {
  PPSTEE_WEIGHTS tmpWgtflag = PPSTEE_WEIGHTS_NONE;

  // set up sets containing non-equipartitioned stage indices
  for (vector<pair<PPSTEE_STAGE, PPSteeWeights*> >::const_iterator it = stages.begin(); it != stages.end(); ++it) {
    if (it->first != PPSTEE_STAGE_REMOVED) {
      if ((it->second->getType() == PPSTEE_WEIGHTS_ALL || it->second->getType() == PPSTEE_WEIGHTS_ONLYVERTEX) && not it->second->areVertexWeightsEquipartitioned()) {
        stagesContainingVertexWeights.insert(it-stages.begin());
      }
      if ((it->second->getType() == PPSTEE_WEIGHTS_ALL || it->second->getType() == PPSTEE_WEIGHTS_ONLYEDGE) && not it->second->areEdgeWeightsEquipartitioned()) {
        stagesContainingEdgeWeights.insert(it-stages.begin());
      }
    }
  }

  // determine available wgtflag
  if (stagesContainingVertexWeights.empty()) {
    if (stagesContainingEdgeWeights.empty()) {
      // all stages are equipartitioned: no need for weights at all
      return PPSTEE_WEIGHTS_NONE;
    } else {
      // only edge weights
      tmpWgtflag = PPSTEE_WEIGHTS_ONLYEDGE;
    }
  } else {
    if (stagesContainingEdgeWeights.empty()) {
      // only vertex weights
      tmpWgtflag = PPSTEE_WEIGHTS_ONLYVERTEX;
    } else {
      // both
      tmpWgtflag = PPSTEE_WEIGHTS_ALL;
    }
  }

  // match available and requested wgtflag (targetWgtflag == PPSTEE_WEIGHTS_NONE was already processed)
  assert(tmpWgtflag != PPSTEE_WEIGHTS_NONE);
  if (tmpWgtflag == PPSTEE_WEIGHTS_ALL) {
    tmpWgtflag = targetWgtflag;
  } else if (targetWgtflag != PPSTEE_WEIGHTS_ALL) { // if it is 'all', we can keep wgtflag as it is.
    if (tmpWgtflag != targetWgtflag) {
      // either 'only vertex vs. only edge' or vice versa, i.e. requested weights are not available (or vice versa): return none
      return PPSTEE_WEIGHTS_NONE;
    } // else (tmpWgtflag == targetWgtflag): everything's fine, keep tmpWgtflag.
  }

  return tmpWgtflag;
}

PPSTEE_WEIGHTS PPStee::getWeightsFlag(PPSTEE_WEIGHTS targetWgtflag) {
  set<int> sv, se;
  return getWeightsFlagAndStageSets(targetWgtflag, sv, se);
}

template<typename T>
PPSTEE_ERROR PPStee::assembleWeights(const int numTargetStages, PPSTEE_WEIGHTS targetWgtflag, const bool contiguousStages, const bool doAlloc, T** vertexWeights, T** edgeWeights, int* numStages, int* wgtflag) {
  // init return values
  if (doAlloc) {
    if (*vertexWeights) { free(*vertexWeights); }
    *vertexWeights = NULL;
    if (*edgeWeights) { free(*edgeWeights); }
    *edgeWeights = NULL;
  }
  *numStages = 0;
  *wgtflag = PPSTEE_WEIGHTS_NONE;

  if (numTargetStages == 0 || getNumberOfActiveStages() == 0 || targetWgtflag == PPSTEE_WEIGHTS_NONE) {
    // numTargetStages == 0: disregard all weights and return empty arrays
    return PPSTEE_ERROR_SUCCESS;
  } else {
    // numTargetStages >= 1:

    // count non-equipartitioned stages and set weight flag
    set<int> stagesContainingVertexWeights;
    set<int> stagesContainingEdgeWeights;
    *wgtflag = getWeightsFlagAndStageSets(targetWgtflag, stagesContainingVertexWeights, stagesContainingEdgeWeights);
    if (*wgtflag == PPSTEE_WEIGHTS_NONE) {
      return PPSTEE_ERROR_SUCCESS;
    }
    set<int> stagesContainingBothWeights(stagesContainingVertexWeights);
    stagesContainingBothWeights.insert(stagesContainingEdgeWeights.begin(), stagesContainingEdgeWeights.end());

    // some stages are not equipartitioned: build vertex and edge weights arrays

    assert(*wgtflag != PPSTEE_WEIGHTS_NONE);
    if (numTargetStages == 1) {
      // numTargetStages == 1: accumulate weights and build single-stage arrays (note: In contrast to PTScotch, Zoltan requires weights data for a specific vertex/edge to be contiguous in memory, however, it does not support multi-phase partitioning (yet))

      if (*wgtflag == PPSTEE_WEIGHTS_ALL || *wgtflag == PPSTEE_WEIGHTS_ONLYVERTEX) {
        gatherWeightsData<T>(PPSTEE_GRAPHDATUM_VERTEX, stagesContainingVertexWeights, true, doAlloc, 0, vertexWeights);
      }

      if (*wgtflag == PPSTEE_WEIGHTS_ALL || *wgtflag == PPSTEE_WEIGHTS_ONLYEDGE) {
        gatherWeightsData<T>(PPSTEE_GRAPHDATUM_EDGE, stagesContainingEdgeWeights, true, doAlloc, 0, edgeWeights);
      }

      if (*wgtflag == PPSTEE_WEIGHTS_ALL) {
        // add vertex weights of stages with edge-only weights. (Currently, just add number. Possible extension: add non-trivial weight value of equipartitioned stage.)
        transform(*vertexWeights, *vertexWeights + graph->getVertloccnt(), *vertexWeights, bind1st(plus<T>(), stagesContainingBothWeights.size() - stagesContainingVertexWeights.size()));

        // add edge weights of stages with vertex-only weights. (Currently, just add number. Possible extension: add non-trivial weight value of equipartitioned stage.)
        transform(*edgeWeights, *edgeWeights + graph->getEdgeloccnt(), *edgeWeights, bind1st(plus<T>(), stagesContainingBothWeights.size() - stagesContainingEdgeWeights.size()));
      }

      // return
      *numStages = 1;
      return PPSTEE_ERROR_SUCCESS;

    } else {
      // numTargetStages >1: build multi-stage arrays

      // calculate #numStages
      *numStages = stagesContainingBothWeights.size();

      if (*wgtflag == PPSTEE_WEIGHTS_ALL || *wgtflag == PPSTEE_WEIGHTS_ONLYVERTEX) {
        gatherWeightsData<T>(PPSTEE_GRAPHDATUM_VERTEX, stagesContainingVertexWeights, false, doAlloc, *numStages - stagesContainingVertexWeights.size(), vertexWeights);
      }

      if (*wgtflag == PPSTEE_WEIGHTS_ALL || *wgtflag == PPSTEE_WEIGHTS_ONLYEDGE) {
        gatherWeightsData<T>(PPSTEE_GRAPHDATUM_EDGE, stagesContainingEdgeWeights, false, doAlloc, *numStages - stagesContainingEdgeWeights.size(), edgeWeights);
      }

      if (*wgtflag == PPSTEE_WEIGHTS_ALL) {
        // append vertex weights of stages with edge-only weights. (Currently, just append number. Possible extension: add non-trivial weight value of equipartitioned stage.)
        fill(*vertexWeights + stagesContainingVertexWeights.size() * graph->getVertloccnt(),
             *vertexWeights + *numStages * graph->getVertloccnt(),
             stagesContainingBothWeights.size() - stagesContainingVertexWeights.size());

        // append edge weights of stages with vertex-only weights. (Currently, just append number. Possible extension: add non-trivial weight value of equipartitioned stage.)
        fill(*edgeWeights + stagesContainingEdgeWeights.size() * graph->getEdgeloccnt(),
             *edgeWeights + *numStages * graph->getEdgeloccnt(),
             stagesContainingBothWeights.size() - stagesContainingEdgeWeights.size());
      }

      if (not contiguousStages) {
        // non-contiguous stages requested: transpose arrays
        if (*wgtflag == PPSTEE_WEIGHTS_ALL || *wgtflag == PPSTEE_WEIGHTS_ONLYVERTEX) {
          transpose<T>(*numStages, graph->getVertloccnt(), vertexWeights);
        }
        if (*wgtflag == PPSTEE_WEIGHTS_ALL || *wgtflag == PPSTEE_WEIGHTS_ONLYEDGE) {
          transpose<T>(*numStages, graph->getEdgeloccnt(), edgeWeights);
        }
      }

      // return
      return PPSTEE_ERROR_SUCCESS;

    }
  }

  // if this line is hit something's terribly wrong!
  return PPSTEE_ERROR_UNKNOWNERROR;
}

PPSTEE_ERROR PPStee::getPartitioningParmetis(PPSteePart** part) {
#ifdef PPSTEE_WITH_PARMETIS
  // grab required data
  int wgtflag; // = getWeightsFlag();
  int numflag = 0;
  int ncon; // = getNumberOfActiveStages();
  int nparts = graph->getMpiN();

  int* vwgt = NULL;
  int* adjwgt = NULL;
  assembleWeights<int>(getNumberOfActiveStages(), PPSTEE_WEIGHTS_ALL, false, true, &vwgt, &adjwgt, &ncon, &wgtflag);
  if (ncon == 0) { ++ncon; } // ParMETIS does need ncon > 0 even if there are no weights to be submitted

  real_t* tpwgts = new real_t[ncon * nparts];
  fill_n(tpwgts, ncon * nparts, 1. / nparts); // for each balance constraint (=stage) every part (=processor) is of equal size

  real_t* ubvec = new real_t[ncon];
  fill_n(ubvec, ncon, 1.05); // imbalance tolerance as recommended in ParMETIS manual

  int* options = new int[1];
  options[0] = 0; // use default values for options

  int edgecut;

  int* newpart = new int[graph->getVertloccnt()];

  MPI_Comm mpi = graph->getMpiComm();

  // call ParMETIS
  ParMETIS_V3_PartKway(graph->getVertglbtab(),
      graph->getVertloctab(),
      graph->getEdgeloctab(),
      vwgt,
      adjwgt,
      &wgtflag,
      &numflag,
      &ncon,
      &nparts,
      tpwgts,
      ubvec,
      options,
      &edgecut,
      newpart,
      &mpi
      );

  // save partitioning
  *part = new PPSteePart(graph, newpart);

  // clean up
  free(vwgt);
  free(adjwgt);
  delete [] tpwgts;
  delete [] ubvec;
  delete [] options;
  delete [] newpart;

  // return
  //return (err == METIS_OK ? PPSTEE_ERROR_SUCCESS : PPSTEE_ERROR_PARMETIS);
  return PPSTEE_ERROR_SUCCESS; // TODO: error handling


#else
  // ParMETIS was not installed with PPStee
  return PPSTEE_ERROR_PARTITIONERNOTINSTALLED;
#endif // PPSTEE_WITH_PARMETIS
}

PPSTEE_ERROR PPStee::getPartitioningPtscotch(PPSteePart** part) {
#ifdef PPSTEE_WITH_PTSCOTCH
  int* vwgt = NULL;
  int* adjwgt = NULL;
  int numStages;
  int wgtFlag;
  assembleWeights<int>(1, PPSTEE_WEIGHTS_ALL, true, true, &vwgt, &adjwgt, &numStages, &wgtFlag);


  SCOTCH_Strat stratptr;
  SCOTCH_stratInit(&stratptr); // default partitioning strategy
/*
  SCOTCH_stratDgraphMapBuild(&stratptr,
      SCOTCH_STRATQUALITY,
      graph->getMpiN(),
      graph->getMpiN(),
      0.01);
*/

  int* newpart = new int[graph->getVertloccnt()];

  // call PTScotch
  SCOTCH_Dgraph* grafptr = SCOTCH_dgraphAlloc();
  SCOTCH_dgraphInit(grafptr, graph->getMpiComm());
  SCOTCH_dgraphBuild(grafptr,
      0,
      graph->getVertloccnt(),
      graph->getVertloccnt(),
      graph->getVertloctab(),
      NULL,
      vwgt,
      NULL,
      graph->getEdgeloccnt(),
      graph->getEdgeloccnt(),
      graph->getEdgeloctab(),
      NULL,
      adjwgt);
  int checkVal = SCOTCH_dgraphCheck(grafptr);
  assert(checkVal==0);

  SCOTCH_dgraphPart(grafptr, graph->getMpiN(), &stratptr, newpart);

  // save partitioning
  *part = new PPSteePart(graph, newpart);

  // clean up
  SCOTCH_dgraphExit(grafptr);
  SCOTCH_dgraphFree(grafptr);
  free(vwgt);
  free(adjwgt);
  delete [] newpart;

  if (checkVal==1) { return PPSTEE_ERROR_PTSCOTCH; } // TODO: error handling
  return PPSTEE_ERROR_SUCCESS; // TODO: error handling

#else
  // PTScotch was not installed with PPStee
  return PPSTEE_ERROR_PARTITIONERNOTINSTALLED;
#endif
}


PPSTEE_ERROR PPStee::getPartitioningZoltan(PPSteePart** part) {
#ifdef PPSTEE_WITH_ZOLTAN
  // init Zoltan
  float version;
  Zoltan_Initialize(NULL, NULL, &version); // argc and argv replaced
  struct Zoltan_Struct* zz;
  zz = Zoltan_Create(graph->getMpiComm());


  // Use only for debug:
  // Zoltan_Memory_Debug(3);
  // Zoltan_Set_Param(zz, "PHG_OUTPUT_LEVEL", "9");
  // Zoltan_Set_Param(zz, "CHECK_HYPERGRAPH", "1");

  // set params
  Zoltan_Set_Param(zz, "LB_METHOD", "GRAPH");
  Zoltan_Set_Param(zz, "LB_APPROACH", "PARTITION"); // "REPARTITION" (and "REFINE") possible
  Zoltan_Set_Param(zz, "RETURN_LISTS", "PARTS");
  // Zoltan_Set_Param(zz, "IMBALANCE_TOL", "1.05"); // default: 1.1
  // NUM_GLOBAL_PARTS defaults to processor count, NUM_LOCAL_PARTS defaults to 1 (number of parts generated)
  // NUM_GID_ENTRIES and NUM_LID_ENTRIES default to 1 (number of ints used to represent Zoltan IDs)

  // set weights
  PPSTEE_WEIGHTS wgtFlag = getWeightsFlag(PPSTEE_WEIGHTS_ALL);
  const char zero[] = "0";
  const char one[] = "1";
  if (wgtFlag == PPSTEE_WEIGHTS_ALL || wgtFlag == PPSTEE_WEIGHTS_ONLYVERTEX) {
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", one);
  } else {
    Zoltan_Set_Param(zz, "OBJ_WEIGHT_DIM", zero);
  }
  if (wgtFlag == PPSTEE_WEIGHTS_ALL || wgtFlag == PPSTEE_WEIGHTS_ONLYEDGE) {
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", one);
  } else {
    Zoltan_Set_Param(zz, "EDGE_WEIGHT_DIM", zero);
  }

  // set callback functions
  Zoltan_Set_Fn(zz, ZOLTAN_NUM_OBJ_FN_TYPE, (void (*)()) &PPStee::zoltanNumObjFn, this); // vertloccnt
  Zoltan_Set_Fn(zz, ZOLTAN_OBJ_LIST_FN_TYPE, (void (*)()) &PPStee::zoltanObjListFn, this); // global_ids[i] = vertglbtab[mpiMe] + i, local_ids[i] = i, i = 0, ..., vertloccnt; obj_wgts: contiguous in vertex weights for each vertex, i.e. obj_wgts = float[vertloccnt][wgt_dim]
  Zoltan_Set_Fn(zz, ZOLTAN_NUM_EDGES_MULTI_FN_TYPE, (void (*)()) &PPStee::zoltanNumEdgesMultiFn, this); // modified vertloctab: num_edges[i] = vertloctab[i+1]-vertloctab[i], i = 0, ..., num_obj-1 (i = local_ids[i] ?)
  Zoltan_Set_Fn(zz, ZOLTAN_EDGE_LIST_MULTI_FN_TYPE, (void (*)()) &PPStee::zoltanEdgeListMultiFn, this); // basically edgeloctab: nbor_global_id = edgeloctab (len = edgeloccnt); nbor_procs[i] = getProcFromVertId(nbor_global_id[i]); ewgts: contiguous in edge weights for each edge, i.e. ewgts = float[edgeloccnt][wgt_dim]

  // prepare partitioning output
  int changes, num_imp, num_exp, *imp_procs, *exp_procs;
  int *imp_to_part, *exp_to_part;
  int num_gid_entries, num_lid_entries;
  ZOLTAN_ID_PTR imp_global_ids, exp_global_ids;
  ZOLTAN_ID_PTR imp_local_ids, exp_local_ids;

  // get partition
  Zoltan_LB_Partition(zz, &changes, &num_gid_entries, &num_lid_entries,
      &num_imp, &imp_global_ids, &imp_local_ids, &imp_procs, &imp_to_part,
      &num_exp, &exp_global_ids, &exp_local_ids, &exp_procs, &exp_to_part); // export_procs = parts

  // save partitioning
  int* newpart = new int[num_exp];
  for (int i=0; i<num_exp; ++i) {
    newpart[i] = exp_procs[i];
  }
  *part = new PPSteePart(graph, newpart);

  // free Zoltan allocated arrays
  Zoltan_LB_Free_Part(&imp_global_ids, &imp_local_ids, &imp_procs, &imp_to_part);
  Zoltan_LB_Free_Part(&exp_global_ids, &exp_local_ids, &exp_procs, &exp_to_part);

  // free Zoltan data structure before ending application.
  Zoltan_Destroy(&zz);

  return PPSTEE_ERROR_SUCCESS; // TODO: error handling

#else
  // Zoltan was not installed with PPStee
  return PPSTEE_ERROR_PARTITIONERNOTINSTALLED;
#endif // PPSTEE_WITH_ZOLTAN
}


int PPStee::zoltanNumObjFn(void* data, int* ierr) {
  PPStee* ppstee = reinterpret_cast<PPStee*>(data); // old c-style cast: PPStee* ppstee = (PPStee*) data;
  *ierr = ZOLTAN_OK;

#ifdef PPSTEE_DEBUG
  std::cout << "[" << ppstee->graph->getMpiMe() << "] " << "ZFn: " << "NumObjFn: " << ppstee->graph->getVertloccnt() << std::endl;
#endif

  return ppstee->graph->getVertloccnt();
}


void PPStee::zoltanObjListFn(void* data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float* obj_wgts, int *ierr) {
  PPStee* ppstee = reinterpret_cast<PPStee*>(data); // old c-style cast: PPStee* ppstee = (PPStee*) data;

#ifdef PPSTEE_DEBUG
  std::cout << "[" << ppstee->graph->getMpiMe() << "] " << "ZFn: " << "ObjListFn: ";
#endif

  for (int i=0; i<ppstee->graph->getVertloccnt(); ++i) {
    global_ids[i] = ppstee->graph->getVertglbtab()[ppstee->graph->getMpiMe()] + i;
    local_ids[i] = i;

#ifdef PPSTEE_DEBUG
    //std::cout << (i==0?"":" ") << ppstee->graph->getVertglbtab()[ppstee->graph->getMpiMe()] + i;
    std::cout << (i==0?"":" ") << global_ids[i];
#endif
  }

#ifdef PPSTEE_DEBUG
  std::cout << std::endl;
#endif

  float* tmpEwgts;
  int numStages, wgtFlag;
  ppstee->assembleWeights<float>(wgt_dim, PPSTEE_WEIGHTS_ONLYVERTEX, false, false, &obj_wgts, &tmpEwgts, &numStages, &wgtFlag);

  *ierr = ZOLTAN_OK;
}


void PPStee::zoltanNumEdgesMultiFn(void *data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges, int *ierr) {
  PPStee* ppstee = reinterpret_cast<PPStee*>(data); // old c-style cast: PPStee* ppstee = (PPStee*) data;

#ifdef PPSTEE_DEBUG
  std::cout << "[" << ppstee->graph->getMpiMe() << "] " << "ZFn: " << "NumEdgesFn: ";
#endif

  for (int i=0; i<num_obj; ++i) {
    num_edges[i] = ppstee->graph->getVertloctab()[i+1] - ppstee->graph->getVertloctab()[i];

#ifdef PPSTEE_DEBUG
    //std::cout << (i==0?"":" ") << ppstee->graph->getVertloctab()[i+1] - ppstee->graph->getVertloctab()[i];
    std::cout << (i==0?"":" ") << num_edges[i];
#endif

  }

#ifdef PPSTEE_DEBUG
  std::cout << std::endl;
#endif

  *ierr = ZOLTAN_OK;
}


void PPStee::zoltanEdgeListMultiFn(void *data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges, ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim, float *ewgts, int *ierr) {
  PPStee* ppstee = reinterpret_cast<PPStee*>(data); // old c-style cast: PPStee* ppstee = (PPStee*) data;

  int offsetGlobalId = ppstee->graph->getVertglbtab()[ppstee->graph->getMpiMe()];
  int sumEdgesSoFar = 0;
  for (int idxObj=0; idxObj < num_obj; ++idxObj) {
    int cur_global_id = global_ids[idxObj];
    int cur_local_id = cur_global_id - offsetGlobalId;
    for (int idxEdge = 0; idxEdge < num_edges[idxObj]; ++idxEdge) {
      nbor_global_id[sumEdgesSoFar + idxEdge] = ppstee->graph->getEdgeloctab()[ppstee->graph->getVertloctab()[cur_local_id] + idxEdge];
    }
    sumEdgesSoFar += num_edges[idxObj];
  }

#ifdef PPSTEE_DEBUG
  std::cout << "[" << ppstee->graph->getMpiMe() << "] " << "ZFn: " << "EdgeListMultiFn: ";
  std::cout << "offGlb = " << offsetGlobalId
    << " num_obj = " << num_obj
    << " sumEdgesSoFar = " << sumEdgesSoFar;
#endif

  int* mappingGlobalIdToProcIdx = (int*) malloc(ppstee->graph->getVertglbtab()[ppstee->graph->getMpiN()] * sizeof(int));
  for (int idxProc = 0; idxProc < ppstee->graph->getMpiN(); ++idxProc) {
    for (int idxGlobalId = ppstee->graph->getVertglbtab()[idxProc]; idxGlobalId < ppstee->graph->getVertglbtab()[idxProc + 1]; ++idxGlobalId) {
      mappingGlobalIdToProcIdx[idxGlobalId] = idxProc;
    }
  }
  for (int idxNbor = 0; idxNbor < sumEdgesSoFar; ++idxNbor) {
    nbor_procs[idxNbor] = mappingGlobalIdToProcIdx[nbor_global_id[idxNbor]];
  }
  free(mappingGlobalIdToProcIdx);

#ifdef PPSTEE_DEBUG
  std::cout << std::endl;
#endif


/*
 *  for (int i=0; i<ppstee->graph->getEdgeloccnt(); ++i) {
 *    nbor_global_id[i] = ppstee->graph->getEdgeloctab()[i];
 *
 *    for (int j=1; j<=ppstee->graph->getMpiN(); ++j) {
 *      if (ppstee->graph->getEdgeloctab()[i] < ppstee->graph->getVertglbtab()[j]) {
 *        nbor_procs[i] = j-1;
 *
 *        break;
 *      }
 *    }
 *  }
 */

  float* tmpObjWgt;
  int numStages, wgtFlag;
  ppstee->assembleWeights<float>(wgt_dim, PPSTEE_WEIGHTS_ONLYEDGE, false, false, &tmpObjWgt, &ewgts, &numStages, &wgtFlag);

  *ierr = ZOLTAN_OK;
}

bool PPStee::isCompatibleGraphAndWeights() {
  for (vector<pair<PPSTEE_STAGE, PPSteeWeights*> >::iterator it = stages.begin(); it != stages.end(); ++it) {
    if (it->first != PPSTEE_STAGE_REMOVED) {
      if (graph->getVertloccnt()!=it->second->getVertloccnt() ||
          graph->getEdgeloccnt()!=it->second->getEdgeloccnt()) {
        return false;
      }
    }
  }

  return true;
}
