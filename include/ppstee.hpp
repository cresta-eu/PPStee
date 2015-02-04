/*! \file ppstee.hpp
 * \brief Main header file.
 *
 * \date 18.09.2012
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

/*! \mainpage Mainpage
 *
 * Welcome to
 * \image html logo.png
 * \image latex logo.png "Logo" width=6cm
 *
 * For further information see \ref documentation (particularly \ref doc_basicusage) or main class PPStee.
 *
 *
 * Comments on the layout are gathered in \ref comments.
 *
 *
 *
 *
 * ToDo:
 * <ul>
 * <li><em>Enable submission of graph data with gaps, i.e. some processors may hold zero vertices. (-> requires new communicator.)</em>
 * <li><em>Additional weight input routine (setWeightsData) for Zoltan query functions.</em>
 * <li><em>PPSteeGraph::checkGraphDataConsistency().</em>
 * <li><em>PPSteeWeights::checkCompatibiliyToGraph().</em>
 * <li><em>Proper error handling.</em>
 * </ul>
 *
 *
 *
 * List of changes:
 * <ul>
 * <li>2013-11-12: <em>v0.3.0:</em>
 *  <ul>
 *  <li><em>Major revision of handling of weights.</em>
 *  <li><em>Added script audit.py for integration of Vera++ and cppcheck.</em>
 *  <li><em>Updated example file.</em>
 *  </ul>
 * <li>2013-08-28: <em>v0.2.1:</em>
 *  <ul>
 *  <li><em>Some minor changes and bugfixes.</em>
 *  </ul>
 * <li>2013-07-29: <em>v0.2.0:</em>
 *  <ul>
 *  <li><em>Added transformation routines for graph data of Zoltan type.</em>
 *  <li><em>Added conversion from local to global graph data.</em>
 *  </ul>
 * <li>2013-06-07:
 *  <ul>
 *  <li><em>Generator for sample graph data added.</em>
 *  <li><em>Added system tests.</em>
 *  <li><em>Added Zoltan partitioning routine.</em>
 *  </ul>
 * <li>2013-04-18: <em>v0.1.1:</em>
 *  <ul>
 *  <li><em>Some minor changes and bugfixes.</em>
 *  </ul>
 * <li>2013-02-28: <em>v0.1.0</em>
 * <li>2012-11-19: <em>Flow chart added. See \ref doc_purpose.</em>
 * <li>2012-09-26: <em>Sections added: \ref comments; this list of changes.</em>
 * </ul>
 */

/*! \page documentation Global documentation
 * \tableofcontents
 * This page gathers arbitrary information. For object details walk through PPStee.
 *
 * \section doc_purpose Purpose of this interface
 *
 * Main purpose of this interface of pre-processing is to optimise overall simulation load-balance.
 * The basic approach is to wrap current partitioning tools (like ParMETIS, PTScotch and Zoltan) in
 * one piece of software to be able to
 * \li use current simulation code with only minor code changes,
 * \li easily test load-balance of the simulation computed by different partitioners (thus test
 * which applied partitioning method fits best for the simulation),
 * \li simplify integration of future partitioning tools,
 * \li incorporate load costs of other simulation parts, like post-processing. (Sure, this could
 * be done in the core calculation code, but modularity has to be favoured in terms of software
 * engineering.)
 *
 * What the interface actually does is
 * \li receive input of graph and weight data,
 * \li call the selected partitioner which computes the partitioning,
 * \li pass this partitioning to the simulation.
 *
 * What the interface does \b not do is
 * \li move \b any data (except for a small amount of private data used for housekeeping).
 *
 * A chart visualising data/information flow:
 * \image html PPStee_flow.png
 * \image latex PPStee_flow.png "flow chart" width=8cm
 *
 * \section doc_basicusage Basic usage
 * Usage will be demonstrated in part_me.cpp. Mainly this is:
 *  \li create and fill PPSteeGraph and several PPSteeWeights
 *  \li create PPStee object and PPStee::submitGraph()
 *  \li computation submits its weights (PPStee::submitNewStageByWeights( * , #PPSTEE_STAGE_COMPUTATION))
 *  \li visualisation submits its weights (PPStee::submitNewStageByWeights( * , #PPSTEE_STAGE_VISUALISATION))
 *  \li PPStee::getPartitioning()
 *
 *  \li PPStee::registerRepartitioningNeed()
 *  \li PPStee::changeStageByWeights(* , #PPSTEE_STAGE_COMPUTATION)
 *  \li PPStee::changeStageByWeights( * , #PPSTEE_STAGE_VISUALISATION)
 *  \li get (re)partitioning (PPStee::getPartitioning())
 *
 * \section doc_graphdata Graph data
 *
 * \subsection doc_graphdata_parmetis ParMETIS graph data layout
 * \code
 * // global (thus duplicated) data
 * int[mpi_n+1] vtxdist;   // global vertex table, = procvrttab
 *
 * // local data
 * int[] xadj;             // local vertex table, = vertloctab
 * int[] adjncy;           // local edge table, = edgeloctab
 * \endcode
 *
 * \remarks
 * Dimensions are not that simple (but obvious):
 * \li dim(xadj) = vtxdist[mpi_me+1]-vtxdist[mpi_me]+1
 * \li dim(adjncy) = xadj[vtxdist[mpi_me+1]-vtxdist[mpi_me]]
 *
 * \subsection doc_graphdata_ptscotch PTScotch graph data layout
 * \code
 * // global (thus duplicated) data
 * int baseval;                         // index base. PPStee consistently uses C++ default (= 0).
 * int vertglbnbr;                      // global vertex count
 * int edgeglbnbr;                      // global edge count
 * int procglbnbr;                      // global processor count, = mpi_n
 * int[procglbnbr+1] procvrttab;        // global vertex table, =vertglbtab, = vtxdist
 * int[procglbnbr] proccnttab;          // table of vertex numbers per processor, = vtxdist[i+1] - vtxdist[i];
 *
 * // local data
 * int vertlocnbr;                      // local vertex count, = vertloccnt
 * int vertlocmax;                      // maximal local vertex number (/index, not count!);
 *                                      // only needed for "graphs without holes in their global numbering", normally = vertlocnbr
 * int[vertlocnbr+1] vertloctab;        // local vertex table
 * int* vendloctab = vertloctab[1];     // local vertex table, points to end indices
 * int[vertlocnbr] vlblloctab = NULL;   // "local vertex label array"; not needed (cf. remark)
 * int edgelocnbr;                      // local edge count; = edgeloccnt; = vertloctab[vertlocnbr];
 * int edgelocsiz = edgelocnbr;         // always true "if the edge array is compact"; (cf. remark)
 * int[edgelocnbr] edgeloctab;          // local edge table
 * int[] edgegsttab = NULL;             // needs analysis of edgeloctab, but not needed (cf. remark)
 * \endcode
 *
 * \remarks
 *  \li \c edgelocsiz: "is lower-bounded by the minimum size of the edge array required to
 *    encompass all used adjacency values; it is therefore at least equal to the maximum of
 *    the \c vendloctab entries, over all local vertices" [source: PTScotch manual]
 *  \li For further information on \c vlblloctab and \c edgegsttab see remark of \ref doc_conversion_from_parmetis_to_ptscotch).
 *  \li Weights (\c veloloctab and \c edloloctab, see \ref doc_weightdata_ptscotch) are also part of the input to PTScotch graph build routine.
 *
 * \subsection doc_conversion_from_parmetis_to_ptscotch Graph data conversion from ParMETIS to PTScotch
 *
 * \code
 * // global (thus duplicated) data
 * baseval    = 0;                     // default C++
 * vertglbnbr = vtxdist[mpi_n];
 * edgeglbnbr = sum( edgelocnbr );     // needs reduce! (cf. remarks)
 * procglbnbr = mpi_n;
 * procvrttab = vtxdist;
 * proccnttab[i] = vtxdist[i+1] - vtxdist[i];
 *
 * // local data
 * vertlocnbr = vtxdist[mpi_me+1] - vtxdist[mpi_me];
 * vertlocmax = vertlocnbr;            // always true for "graphs without holes in their global numbering"
 * vertloctab = xadj;                  // both of dimension (vertlocnbr + 1)
 * vendloctab = vertloctab + 1;        // implicit dimension equals vertlocnbr
 * veloloctab = vwgt;                  // dimension = vertlocnbr
 * vlblloctab = NULL;                  // "local vertex label array"; dimension = vertlocnbr; not needed (cf. remark)
 * edgelocnbr = xadj[vertlocnbr];
 * edgelocsiz = edgelocnbr;            // always true "if the edge array is compact"
 * edgeloctab = adjncy;                // dimension = edgelocnbr
 * edgegsttab = NULL;                  // needs analysis of edgeloctab, but not needed (cf. remark)
 * edloloctab = adjwgt;                // dimension = edgelocnbr
 * \endcode
 *
 * \remarks
 * \li In PTScotch graph build routine \c vendloctab, \c veloloctab, \c vlblloctab,
 * \c edloloctab and \c edgegsttab are optional. A NULL pointer can be provided. Global data (except for \c baseval) is not needed at all.
 * \li \c edgeglbnbr needs global communication. (However, keep in mind that it is not needed for graph building.)
 * \li MPI is needed for \c vertglbnbr, \c procglbnbr and \c vertlocnbr.
 *
 * \subsection doc_conversion_from_ptscotch_to_parmetis Graph data conversion from PTScotch to ParMETIS
 * \code
 * // global (thus duplicated) data
 * vtxdist = procvrttab;
 *
 * // local data
 * xadj = vertloctab;
 * adjncy = edgeloctab;
 * \endcode
 *
 * \subsection doc_conversion_from_parmetis_like_to_zoltan Graph data conversion from ParMETIS-like data to Zoltan query functions
 * \li ZOLTAN_NUM_OBJ_FN_TYPE: return vertloccnt
 * \li ZOLTAN_OBJ_LIST_FN_TYPE: global_ids[i] = vertglbtab[mpiMe] + i, local_ids[i] = i, i = 0, ..., vertloccnt;\n obj_wgts: contiguous in vertex weights for each vertex, i.e. obj_wgts = float[vertloccnt][wgt_dim]
 * \li ZOLTAN_NUM_EDGES_MULTI_FN_TYPE: modified vertloctab: num_edges[i] = vertloctab[i+1]-vertloctab[i], i = 0, ..., num_obj-1 (i = local_ids[i])
 * \li ZOLTAN_EDGE_LIST_MULTI_FN_TYPE: basically edgeloctab: nbor_global_id = edgeloctab (len = edgeloccnt);\n nbor_procs[i] = getProcFromVertId(nbor_global_id[i]);\n ewgts: contiguous in edge weights for each edge, i.e. ewgts = float[edgeloccnt][wgt_dim]
 *
 * \subsection doc_graphdata_conciliation Graph data conciliation: unified graph data layout
 * Compare private data fields of PPSteeGraph.
 * \code
 * // MPI related
 * int mpi_me;                  // MPI_Comm_rank
 * int mpi_n;                   // MPI_Comm_size
 *
 * // global (thus duplicated) data
 * int[mpi_n+1]    vertglbtab;  // global vertex table; = vtxdist; = procvrttab
 *
 * // local data
 * int               vertloccnt;  // local vertex count
 * int               edgeloccnt;  // local edge count
 * int[vertloccnt+1] vertloctab;  // local vertex table; = xadj
 * int[edgeloccnt]   edgeloctab;  // local edge table; = adjncy
 * \endcode
 *
 * Properties:
 * \li Almost minimal.\n Memory overhead: 2 \c ints (used often). (cf. \ref doc_graphdata_parmetis)
 * \li Little calculation needed for conversions. (cf. \ref doc_conversion_from_parmetis_to_ptscotch and \ref doc_conversion_from_ptscotch_to_parmetis)
 * \li Only global edge count is not stored. Can be computed by reducing local edge counts. However, it is rarely needed.
 *
 *
 * \section doc_weightdata Weight data
 * Compare private data fields of PPSteeWeights.
 *
 * Keep in mind that weights' order is bound to vertices' order.
 *
 * \subsection doc_weightdata_parmetis ParMETIS weight data layout
 * \code
 * // local data
 * int[vertloccnt] vwgt;                // vertex weights
 * int[edgeloccnt] adjwgt;              // edge weights
 * \endcode
 *
 * \remarks
 * Dimensions in ParMETIS' variables:
 * \li dim(vwgt) = vtxdist[mpi_me+1]-vtxdist[mpi_me] = dim(xadj)-1
 * \li dim(adjwgt) = xadj[vtxdist[mpi_me+1]-vtxdist[mpi_me]] = dim(adjncy)
 *
 *
 * \subsection doc_weightdata_ptscotch PTScotch weight data layout
 * \code
 * // local data
 * int[vertlocnbr] veloloctab;          // vertex weights
 * int[edgelocnbr] edloloctab;          // edge weights
 * \endcode
 *
 * \subsection doc_weightdata_zoltan Zoltan weight data
 * Zoltan's weight data is integrated in query functions ZOLTAN_OBJ_LIST_FN_TYPE and ZOLTAN_EDGE_LIST_MULTI_FN_TYPE
 * (cf. \ref doc_conversion_from_parmetis_like_to_zoltan).\n
 * \li ZOLTAN_OBJ_LIST_FN_TYPE: obj_wgts: contiguous in vertex weights for each vertex, i.e. obj_wgts = float[vertloccnt][wgt_dim]
 * \li ZOLTAN_EDGE_LIST_MULTI_FN_TYPE: ewgts: contiguous in edge weights for each edge, i.e. ewgts = float[edgeloccnt][wgt_dim]
 *
 * \subsection doc_weightdata_conciliation Weight data conciliation: unified weight data layout
 * Compare private data fields of PPSteeWeights.
 * \code
 * // local data
 * int[vertloccnt] vertwgttab;          // vertex weights; = vwgt; = veloloctab
 * int[edgeloccnt] edgewgttab;          // edge weights; = adjwgt; = edloloctab
 * \endcode
 *
 * Properties:
 * \li \b Minimal!
 * \li \b No calculation needed for conversions. (cf. \ref doc_weightdata_parmetis and \ref doc_weightdata_ptscotch)
 *
 */

/*! \page comments Questions and answers
 * \tableofcontents
 *
 * \section comments_christoph By Christoph
 *
 * \subsection comments_christoph_0 [Q0.0] Get I it right that the interface intends to read in and write out input data (Matrices) of different formats as well as allow manipulations on the CRESTA internal format?
 * Until now, simulations are responsible for reading graph/matrix data. Thus the interface will not
 * interfere with data in any way. See also \ref comments_christoph_2 and \ref comments_christoph_3.
 *
 * However, some kind of reading mechanism could be a future feature.
 *
 * \subsection comments_christoph_1 [Q0.1] C++? May we need C or Fortran support for this interface, too?
 * In my opinion at least general feasibility is established best in C++ (due to readibility,
 * architecture design and many more). Once this is done, we can think of how much of it should be
 * done in C.
 *
 * Some other points:
 * \li Most partitioners provide C libs. Therefore main interface logic could be done in
 * core functions implemented in C. Or:
 * \li A wrapper to C is not simple, but doable (even afterwards).
 * \li Our main target applications (HemeLB, OpenFOAM) are C++.
 *
 * \subsection comments_christoph_2 [Q0.2] What about data types (e.g. structs) instead of all these global variables? This could help if different data sets have to be mixed in pre/postrocessing.
 * See \ref doc_purpose. Global variables are only pointers to already present data. And: they are
 * determined by the various partitioning tools used.
 *
 * \subsection comments_christoph_3 [Q0.3] What are the actual operations which will be performed on (all) the data?
 * None. See \ref doc_purpose.
 *
 */

#ifndef PPSTEE_HPP
#define PPSTEE_HPP

#include <mpi.h>


#ifdef PPSTEE_WITH_PARMETIS
#include <metis.h>
#include <parmetis.h>
#endif // PPSTEE_WITH_PARMETIS

#ifdef PPSTEE_WITH_PTSCOTCH
#include <ptscotch.h>
#endif // PPSTEE_WITH_PTSCOTCH

#ifdef PPSTEE_WITH_ZOLTAN
#include <zoltan.h>
#endif // PPSTEE_WITH_ZOLTAN



#ifdef PPSTEE_DEBUG
#include <iostream>
#endif



#define PPSTEE_TRANSPOSE_BLOCKSIZE 16  //!< Blocksize for transpose.

#define PPSTEE_INT_SIZE 4  //!< PPStee's internal int size. TODO: Use partitioner specific int sizes.


#include<stdlib.h>
#include<cstring>
#include<vector>
#include<set>
#include<algorithm>
#include<functional>
#include<assert.h>

using std::vector;
using std::set;
using std::pair;
using std::fill_n;
using std::copy;
using std::memcpy;
using std::transform;
using std::plus;
using std::bind1st;
using std::fill;


//// forward declarations
class PPSteeGraph;
class PPSteeWeights;
class PPSteePart;
class PPSteeError;


#include "ppstee_objects.hpp"
#include "ppstee_error.hpp"


//! PPStee stage type enumeration.
/*!
 * Simulation stages can/should pass their extent of workload and communication to pre-processing.
 * These stages are organised in types defining their position in the computation sequence.
 *
 * Submitted weights are grouped by type and within each group arranged by time of submission
 * (first in, first out).
 *
 * The chronological order of types is set to: computation > visualisation > other.
 */
typedef enum ENUM_PPSTEE_STAGETYPE {
  PPSTEE_STAGE_COMPUTATION,    //!< PPSTEE_STAGE_COMPUTATION
  PPSTEE_STAGE_VISUALISATION,  //!< PPSTEE_STAGE_VISUALISATION
  PPSTEE_STAGE_OTHER,          //!< PPSTEE_STAGE_OTHER
  PPSTEE_STAGE_REMOVED         //!< PPSTEE_STAGE_REMOVED
} PPSTEE_STAGE;


//! PPStee base class
/*!
 * This class provides an interface to steer pre-processing.
 *
 * Supported partitioning libraries are ParMETIS and PTScotch. (later: Zoltan)
 *
 * For information on basic usage see \ref doc_basicusage.
 */
class PPStee {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Constructor.
  /*!
   * Private #useOldPartitioning defaults to false.
   */
  PPStee();


  ///@}
  //! \name Submit graph
  ///@{

  //! Submit graph.
  /*!
   * \note If a graph was already submitted isCompatibleGraphAndWeights() is called.
   *
   * \param graph
   *    Graph to be used for partitioning.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR submitGraph(PPSteeGraph& graph);


  ///@}
  //! \name Submit new weights
  ///@{

  //! Submit new stage; redirects to submitNewStageByWeights().
  /*!
   *
   * \param[in] weights
   *    Weights for this new stage. Must match supplied graph.
   * \param[in] type
   *    Stage type of this new stage.
   * \param[out] stageNumber
   *    Number of stage submitted.
   *    (optional, defaults to NULL)
   * \return
   *    Error code.
   */
  PPSTEE_ERROR submitNewStage(PPSteeWeights& weights, PPSTEE_STAGE type, int* stageNumber = NULL) { return submitNewStageByWeights(weights, type, stageNumber); };

  //! Submit new stage.
  /*!
   *
   * \param[in] weights
   *    Weights for this new stage. Must match supplied graph (use default PPSteeWeights constructor)
   * \param[in] type
   *    Stage type of this new stage.
   * \param[out] stageNumber
   *    Number of stage submitted.
   *    (optional, defaults to NULL)
   * \return
   *    Error code.
   */
  PPSTEE_ERROR submitNewStageByWeights(PPSteeWeights& weights, PPSTEE_STAGE type, int* stageNumber = NULL);


  ///@}
  //! \name Change existing weights
  ///@{

  //! Change existing stage; redirects to changeStageByWeights().
  /*!
   *
   * \param[in] weights
   *    Weights for existing stage. Must match supplied graph.
   * \param[in] stageNumber
   *    Number of stage to change.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR changeStage(PPSteeWeights& weights, int stageNumber) { return changeStageByWeights(weights, stageNumber); };

  //! Change existing stage.
  /*!
   *
   * \param[in] weights
   *    Weights for existing stage. Must match supplied graph.
   * \param[in] stageNumber
   *    Number of stage to change.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR changeStageByWeights(PPSteeWeights& weights, int stageNumber);


  ///@}
  //! \name Remove existing weights
  ///@{

  //! Remove existing stage.
  /*!
   * Removing a stage does not alter numbering of stages.
   *
   * \param[in] stageNumber
   *    Stage number to be removed.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR removeStage(int stageNumber);


  ///@}
  //! \name Partitioning
  ///@{


  //! Compute partitioning.
  /*!
   *
   * \param[out] part
   *    Stores computed partitioning.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR getPartitioning(PPSteePart** part);


  ///@}
  /*!
   * \name Repartitioning
   *
   * A (re-)partitioning is faster if it is not computed from scratch but rather built up on an old one.
   * So there are two complementary possibilities:
   *
   * - registerRepartitioningNeed() sets the #useOldPartitioning flag: the next call to getPartitioning()
   * will save the computed partitioning and will make use of it in the following (i.e., the second) call to getPartitioning().
   *
   * - Also, an old partitioning can be submitted beforehand (cf. submitOldPartitioning())
   * resulting in usage of this partitioning already in the first following getPartitioing().
   *
   * Nevertheless, new stages can be submitted (cf. submitNewStage())
   * and existing ones can be changed (cf. changeStage()).
   */
  ///@{

  //! Register need for repartitioning. Just marks the next partitioning to be saved.
  void registerRepartitioningNeed() {useOldPartitioning = true; };

  //! Submit old partitioning to be used for next getPartitioning().
  PPSTEE_ERROR submitOldPartitioning(PPSteePart& part);

  //! Reset need for repartitioning. Next partitioning will not be saved.
  /*!
   * Resets #useOldPartitioning so that the next call to getPartitioning() will not make use of the old partitioning.
   *
   * \note Does \em not touch a possibly present partitioning.
   */
  void unregisterRepatitioningNeed() {useOldPartitioning = false; };


  ///@}
  //! \name Getters/Setters
  ///@{

  //! Get number of active (i.e., not removed) stages.
  int getNumberOfActiveStages() const;

  //! Get number of all (i.e., including removed) stages.
  int getNumberOfAllStages() const {return stages.size(); };

  //! Get stage type.
  PPSTEE_STAGE getStageType(int stage) const {return stages[stage].first; };

  //! Get #oldPart.
  PPSteePart* getOldPart() const {return oldPart; };

  //! Get #useOldPartitioning.
  bool getUseOldPartitioning() const {return useOldPartitioning; };

  ///@}

private:

  //! \name private data fields
  ///@{

  PPSteeGraph* graph; //!< Submitted graph.

  vector<pair<PPSTEE_STAGE, PPSteeWeights*> > stages; //!< Array of submitted stages and their weights.

  PPSteePart* oldPart; //!< Submitted or computed old partitioning.

  bool useOldPartitioning /*= false*/;   //!< Switch to use old partitioning. Useful if repartitioning is necessary.


  ///@}
  //! \name Partitioning
  ///@{


  //! Compute partitioning with ParMETIS
  /*!
   *
   * \param[out] part
   *    Stores computed partitioning.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR getPartitioningParmetis(PPSteePart** part);


  //! Compute partitioning with PTScotch
  /*!
   *
   * \param[out] part
   *    Stores computed partitioning.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR getPartitioningPtscotch(PPSteePart** part);


  //! Compute partitioning with Zoltan
  /*!
   *
   * \param[out] part
   *    Stores computed partitioning.
   * \return
   *    Error code.
   */
  PPSTEE_ERROR getPartitioningZoltan(PPSteePart** part);


  ///@}
  //! \name Zoltan callback functions
  ///@{

  //! Provide number of objects currently assigned to the processor.
  /*!
   * \param[in] data Pointer to user-defined data.
   * \param[out] ierr Error code to be set by function.
   *
   * \return number of objects currently assigned to the processor.
   *
   */
  static int zoltanNumObjFn(void* data, int* ierr);


  //! Provide global ids, local ids and weights of objects currently assigned to the processor.
  /*!
   * \param[in] data Pointer to user-defined data.
   * \param[in] num_gid_entries The number of array entries used to describe a single global ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES.
   * \param[in] num_lid_entries The number of array entries used to describe a single local ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES. (It should be zero if local ids are not used.)
   * \param[out] global_ids Upon return, an array of unique global IDs for all objects assigned to the processor.
   * \param[out] local_ids Upon return, an array of local IDs, the meaning of which can be determined by the application, for all objects assigned to the processor. (Optional.)
   * \param[in] wgt_dim The number of weights associated with an object (typically 1), or 0 if weights are not requested. This value is set through the parameter OBJ_WEIGHT_DIM.
   * \param[out] obj_wgts Upon return, an array of object weights. Weights for object i are stored in obj_wgts[(i-1)*wgt_dim:i*wgt_dim-1].  If wgt_dim=0, the return value of obj_wgts is undefined and may be NULL.
   * \param[out] ierr Error code to be set by function.
   *
   */
  static void zoltanObjListFn(void* data, int num_gid_entries, int num_lid_entries, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int wgt_dim, float* obj_wgts, int *ierr);


  //! Provide number of edges for each object currently assigned to the processor.
  /*!
   * \param[in] data Pointer to user-defined data.
   * \param[in] num_gid_entries The number of array entries used to describe a single global ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES.
   * \param[in] num_lid_entries The number of array entries used to describe a single local ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES. (It should be zero if local ids are not used.)
   * \param[in] num_obj The number of object IDs in arrays global_ids and local_ids.
   * \param[in] global_ids Array of global IDs of objects whose number of edges should be returned.
   * \param[in] local_ids Array of local IDs of objects whose number of edges should be returned. (Optional.)
   * \param[out] num_edges Upon return, an array containing numbers of edges. For object i (specified by global_ids[i*num_gid_entries] and local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), the number of edges should be stored in num_edges[i].
   * \param[out] ierr Error code to be set by function.
   *
   */
  static void zoltanNumEdgesMultiFn(void *data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges, int *ierr);


  //! Provide edge data (ids of neighbored objects and the processor id where they reside) for each object currently assigned to the processor.
  /*!
   * \param[in] data Pointer to user-defined data.
   * \param[in] num_gid_entries The number of array entries used to describe a single global ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES.
   * \param[in] num_lid_entries The number of array entries used to describe a single local ID.  This value is the maximum value over all processors of the parameter NUM_GID_ENTRIES. (It should be zero if local ids are not used.)
   * \param[in] num_obj The number of object IDs in arrays global_ids and local_ids.
   * \param[in] global_ids Array of global IDs of objects whose edge lists should be returned.
   * \param[in] local_ids Array of local IDs of objects whose edge lists should be returned. (Optional.)
   * \param[in] num_edges An array containing numbers of edges for each object in global_ids. For object i (specified by global_ids[i*num_gid_entries] and local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), the number of edges is stored in num_edges[i].
   * \param[out] nbor_global_id Upon return, an array of global IDs of objects sharing edges with the objects specified in global_ids. For object i (specified by global_ids[i*num_gid_entries] and local_ids[i*num_lid_entries], i=0,1,...,num_obj-1), edges are stored in nbor_global_id[sum*num_gid_entries] to nbor_global_id[(sum+num_edges[i])*num_gid_entries-1], where sum = the sum of num_edges[j] for j=0,1,...,i-1.
   * \param[out] nbor_procs Upon return, an array of processor IDs that identifies where the neighboring objects reside. For neighboring object i (stored in nbor_global_id[i*num_gid_entries]), the processor owning the neighbor is stored in nbor_procs[i].
   * \param[in] wgt_dim The number of weights associated with an edge (typically 1), or 0 if edge weights are not requested. This value is set through the parameter EDGE_WEIGHT_DIM.
   * \param[out] ewgts Upon return, an array of edge weights, where ewgts[i*wgt_dim:(i+1)*wgt_dim-1]
   * corresponds to the weights for the ith edge. If wgt_dim=0, the return value of ewgts is undefined and may be NULL.
   * \param[out] ierr Error code to be set by function.
   *
   */
  static void zoltanEdgeListMultiFn(void *data, int num_gid_entries, int num_lid_entries, int num_obj, ZOLTAN_ID_PTR global_ids, ZOLTAN_ID_PTR local_ids, int *num_edges, ZOLTAN_ID_PTR nbor_global_id, int *nbor_procs, int wgt_dim, float *ewgts, int *ierr);


  ///@}
  //! \name Helpers
  ///@{

  //! Checks compatibility of graph and weights.
  bool isCompatibleGraphAndWeights();

  /**
   * \brief Transpose.
   *
   * Transpose array using a blocking scheme (cf. #PPSTEE_TRANSPOSE_BLOCKSIZE).
   * \param[in] h Height.
   * \param[in] w Width.
   * \param[in,out] array Array.
   *
   * @return Error code.
   */
  template <typename T>
  PPSTEE_ERROR transpose(const int h, const int w, T** array);

  /**
   * \brief Determines weight flag based on target weight flag and fills sets with indices of stages containing non-equipartitioned weights.
   *
   * \param[in] targetWgtflag Specifies requested type of weights. (cf. PPStee::assembleWeights() for more details.)
   * \param[in] stagesContainingVertexWeights Set of stage indices containing non-equipartitioned vertex weights.
   * \param[in] stagesContainingEdgeWeights Set of stage indices containing non-equipartitioned edge weights.
   *
   * @return Resulting weight flag.
   */
  PPSTEE_WEIGHTS getWeightsFlagAndStageSets(PPSTEE_WEIGHTS targetWgtflag, set<int>& stagesContainingVertexWeights, set<int>& stagesContainingEdgeWeights);

  /**
   * \brief Determines weight flag based on target weight flag.
   *
   * \param[in] targetWgtflag Specifies requested type of weights. (cf. PPStee::assembleWeights() for more details.)
   *
   * @return Resulting weight flag.
   */
  PPSTEE_WEIGHTS getWeightsFlag(PPSTEE_WEIGHTS targetWgtflag);

  /**
   * \brief Gathers weights data.
   *
   * \param[in] type Gather weights of this type, i.e., either PPSTEE_GRAPHDATUM_VERTEX or PPSTEE_GRAPHDATUM_EDGE.
   * \param[in] lstStages Gather weights of these stages.
   * \param[in] accumulate Accumulate weights of provided stages in a single stage (true) or build multi-stage array (false).
   * \param[in] doAlloc Specifies if array should be allocated.
   * \param[in] numSpareStages If an multi-stage array is requested, alloc #\paramname{numSpareStages} stages additionally. (Do \em not forget to set \paramname{doAlloc}.)
   * \param[out] array Resulting array.
   *
   * @return Error code.
   */
  template <typename T>
  PPSTEE_ERROR gatherWeightsData(PPSTEE_GRAPHDATUM type, const set<int>& lstStages, const bool accumulate, const bool doAlloc, const int numSpareStages, T** array);

  /**
   * \brief Assembles weights arrays.
   *
   * \param[in] numTargetStages Specifies resulting number of stages. For details, see below.
   * \param[in] targetWgtflag Specificies requested type of weights. This marks the maximum-possible weights extent and therefore does \em not have to be the resulting type of weights. If, e.g., only vertex weights are requested, but there are no vertex weights available (i.e., only edge weights are provided), the resulting \paramname{wgtflag} will be none.
   * \param[in] contiguousStages Specifies which part of data should be aligned: either stage weights are contiguous or each vertex/edge weights (the former is required by PTScotch, the latter by Zoltan and ParMETIS).
   * \param[in] doAlloc Specifies if arrays should be allocated.
   * \param[out] vertexWeights Vertex weights.
   * \param[out] edgeWeights Edge weights.
   * \param[out] numStages Resulting number of stages.
   * \param[out] wgtflag wgtflag as required (and defined) by ParMETIS.
   *
   * Options for \paramname{numTargetStages}:
   * - \paramname{numTargetStages} = 0: Disregards all weights and returns empty arrays; this is needed for Zoltan (cf. zoltanObjListFn() or zoltanEdgeListMultiFn()))
   * - \paramname{numTargetStages} = 1: Accumulates all stages in a single stage.
   * - \paramname{numTargetStages} > 1: Counts and uses all stages provided and returns #\paramname{numStages} stages, disregarding the specific value of \paramname{numTargetStages}.
   *
   * \par
   * \note Regardless of the specific value of \paramname{numTargetStages}, numStages \em may differ.
   *
   * \par
   * \note In contrast to PTScotch, Zoltan and ParMETIS require weights data for a specific vertex/edge to be contiguous in memory.
   *
   * \par
   * \note PTScotch and Zoltan do not support multi-phase partitioning (yet).
   *
   * @return Error code.
   */
  template <typename T>
  PPSTEE_ERROR assembleWeights(const int numTargetStages, PPSTEE_WEIGHTS targetWgtflag, const bool contiguousStages, const bool doAlloc, T** vertexWeights, T** edgeWeights, int* numStages, int* wgtflag);


  ///@}
};


#endif /* PPSTEE_HPP_ */
