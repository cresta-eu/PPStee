/*! \file ppstee_objects.hpp
 * \brief Defining objects used by PPStee.
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

#ifndef PPSTEE_OBJECTS_HPP
#define PPSTEE_OBJECTS_HPP


//! PPStee graph type enumeration.
/*!
 * Several (derived) graph types are available. Their types are defined here.
 */
typedef enum ENUM_PPSTEE_GRAPHTYPE {
  PPSTEE_GRAPH_UNDETERMINED,  //!< PPSTEE_GRAPH_UNDETERMINED
  PPSTEE_GRAPH_PARMETIS,      //!< PPSTEE_GRAPH_PARMETIS
  PPSTEE_GRAPH_PTSCOTCH,      //!< PPSTEE_GRAPH_PTSCOTCH
  PPSTEE_GRAPH_ZOLTAN         //!< PPSTEE_GRAPH_ZOLTAN
} PPSTEE_GRAPH;

//! PPStee weights type enumeration.
/*!
 * Weights type specifies usage of vertex, edge or both weights.
 */
typedef enum ENUM_PPSTEE_WEIGHTSTYPE {
  PPSTEE_WEIGHTS_NONE = 0, //!< PPSTEE_WEIGHTS_NONE
  PPSTEE_WEIGHTS_ALL = 3, //!< PPSTEE_WEIGHTS_ALL (default)
  PPSTEE_WEIGHTS_ONLYVERTEX = 2,  //!< PPSTEE_WEIGHTS_ONLYVERTEX
  PPSTEE_WEIGHTS_ONLYEDGE = 1 //!< PPSTEE_WEIGHTS_ONLYEDGE
} PPSTEE_WEIGHTS;

//! PPStee data access type enumeration.
/*!
 * Specifies the type of data access PPStee performs.
 *
 * \note This refers only to arrays used, i.e.: {vert,edge}{glb,loc,wgt}tab.
 */
typedef enum ENUM_PPSTEE_DATAACCESSTYPE {
  PPSTEE_DATAACCESS_COPY, //!< PPSTEE_DATAACCESS_COPY
  PPSTEE_DATAACCESS_VIEW, //!< PPSTEE_DATAACCESS_VIEW
  PPSTEE_DATAACCESS_QRYFN //!< PPSTEE_DATAACCESS_QRYFN
} PPSTEE_DATAACCESS;

//! PPStee graph datum type enumeration.
/*!
 * Specifies the type of datum in a graph.
 */
typedef enum ENUM_PPSTEE_GRAPHDATUMTYPE {
  PPSTEE_GRAPHDATUM_VERTEX, //!< PPSTEE_GRAPHDATUM_VERTEX
  PPSTEE_GRAPHDATUM_EDGE //!< PPSTEE_GRAPHDATUM_EDGE
} PPSTEE_GRAPHDATUM;


#include<stdlib.h>
#include<cstring> // for memcpy


#include "ppstee.hpp"
#include "ppstee_error.hpp"

using std::memcpy;
using std::fill;


//! PPStee graph class.
/*!
 * Graphs are basic information objects to supply to partitioning class PPStee.
 */
class PPSteeGraph {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Constructor.
  PPSteeGraph(MPI_Comm mpiComm, PPSTEE_DATAACCESS access = PPSTEE_DATAACCESS_COPY);

  //! Copy/convert constructor.
  PPSteeGraph(PPSteeGraph* existingGraph, PPSTEE_GRAPH targetType);

  //! Destructor.
  ~PPSteeGraph();


  ///@}
  //! \name Mesh changes
  ///@{

  //! Refine underlying mesh automatically.
  /*!
   * An automated mesh refinement could be implemented here.
   *
   * \warning Not implemented yet.
   *
   * \return
   *    Error code.
   */
  PPSTEE_ERROR refineMeshAutomatically();


  ///@}
  //! \name Consistency checks
  ///@{

  //! Check consistency of graph data.
  /*!
   * Analyse graph data and report when inconsistent.
   */
  PPSTEE_ERROR checkGraphDataConsistency();


  ///@}
  //! \name Free data
  ///@{

  //! Frees memory of private arrays holding most of this class' data.
  /*!
   * \warning Use with care.
   */
  void freeMemoryIfHoldingDataCopies();


  ///@}
  //! \name Getters/Setters
  ///@{

  //! Get #type.
  PPSTEE_GRAPH getType() const {return type; };

  //! Get #dataAccess.
  PPSTEE_DATAACCESS getDataAccess() const {return dataAccess; };

  //! Get #onlyLocalData.
  bool getOnlyLocalData() const {return onlyLocalData; };

  //! Get #mpiComm.
  MPI_Comm getMpiComm() const {return mpiComm; };
  //! Get #mpiMe.
  int getMpiMe() const {return mpiMe; };
  //! Get #mpiN.
  int getMpiN() const {return mpiN; };

  //! Get #vertglbtab.
  virtual int* getVertglbtab();

  //! Get #vertloccnt.
  int getVertloccnt() const {return vertloccnt; };
  //! Get #edgeloccnt.
  int getEdgeloccnt() const {return edgeloccnt; };
  //! Get #vertloctab.
  virtual int* getVertloctab() const {return vertloctab; };
  //! Get #edgeloctab.
  virtual int* getEdgeloctab() const {return edgeloctab; };

  ///@}

protected:

  //! \name Graph data and property access
  ///@{

  PPSTEE_ERROR setGraphData(int* vertglbtab, int* vertloctab, int* edgeloctab);  //!< Set graph data. \note Takes care of #dataAccess.
  PPSTEE_ERROR setLocalGraphData(int vertloccnt, int* vertloctab, int* edgeloctab);  //!< Set only local graph data; specifically, #vertglbtab is not set nor used. See PPStee::getPartitioningPtscotch(). \note Takes care of #dataAccess.
  void setType(PPSTEE_GRAPH type);  //!< Set #type of graph.
  PPSTEE_ERROR buildGlobalDataFromLocalData();  //!< Uses #vertloccnt and MPI_Allgather to assemble #vertglbtab.


  ///@}

protected:

  //! \name Protected data fields
  ///@{

  PPSTEE_GRAPH type;  //!< Type of this graph.

  PPSTEE_DATAACCESS dataAccess; //!< Type of data access for this graph.

  bool onlyLocalData; //!< Specifies if only local data is present, i.e. #vertglbtab is not set (yet).

  MPI_Comm mpiComm; //!< MPI_Comm
  int mpiMe;  //!< MPI_Comm_rank
  int mpiN; //!< MPI_Comm_size

  int* vertglbtab;  //!< global vertex table; = vtxdist; = procvrttab

  int vertloccnt; //!< local vertex count
  int edgeloccnt;  //!< local edge count
  int* vertloctab;  //!< local vertex table; = xadj
  int* edgeloctab;  //!< local edge table; = adjncy


  ///@}

};


//! PPStee ParMETIS graph class
/*!
 * Graph class for access in native ParMETIS description.
 */
class PPSteeGraphParmetis : public PPSteeGraph {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Constructor. Calls create().
  PPSteeGraphParmetis(MPI_Comm mpiComm, int* vtxdist, int* xadj, int* adjncy, PPSTEE_DATAACCESS access = PPSTEE_DATAACCESS_VIEW);

  //! Convert constructor.
  /*!
   * Construct ParMETIS graph based on another graph. Especially, this includes graphs of other types.
   *
   * \param[in] existingGraph
   */
  PPSteeGraphParmetis(PPSteeGraph* existingGraph);


  ///@}
  //! \name Graph operations.
  ///@{

  //! Create ParMETIS graph.
  /*!
   *
   * \param[in] vtxdist,xadj,adjncy arrays needed for a ParMETIS graph
   * \return
   *    Error code.
   */
  PPSTEE_ERROR create(int* vtxdist, int* xadj, int* adjncy);

  ///@}
  //! \name Getters/Setters
  ///@{

  ///@}

//private:

  //! \name Private data fields
  ///@{

  ///@}
};

//! PPStee PTScotch graph class
/*!
 * Graph class for access in native PTScotch description.
 */
class PPSteeGraphPtscotch : public PPSteeGraph {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Full-parameter constructor. Calls create().
  /*!
   * All parameters are identical in type and order to a call to SCOTCH_dgraphBuild().
   */
  PPSteeGraphPtscotch(MPI_Comm mpiComm,
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
                      PPSTEE_DATAACCESS access = PPSTEE_DATAACCESS_COPY
                      );

  //! Minimal-parameter constructor. Calls create().
  /*!
   * Only needed parameters has to be provided; redundant parameters of the full list (cf. PPSteeGraphPtscotch(int, int, int, int*, int*, int*, int*, int, int, int*, int*, int*)) were stripped off. Names match those in SCOTCH_dgraphBuild().
   */
  PPSteeGraphPtscotch(MPI_Comm mpiComm,
                      int vertlocnbr,
                      int* vertloctab,
                      int* edgeloctab,
                      PPSTEE_DATAACCESS access = PPSTEE_DATAACCESS_COPY
                      );


  //! Convert constructor.
  /*!
   * Construct PTScotch graph based on another graph. Especially, this includes graphs of other types.
   *
   * Conversion from ParMETIS to PTScotch: see \ref doc_conversion_from_parmetis_to_ptscotch
   *
   * \param[in] existingGraph
   */
  PPSteeGraphPtscotch(PPSteeGraph* existingGraph);


  ///@}
  //! \name Graph operations.
  ///@{

  //! Create PTScotch graph.
  /*!
   *
   * \param[in] vertlocnbr, vertloctab, edgeloctab Arrays (really) needed for a PTScotch graph. (cf. \ref doc_graphdata_ptscotch)
   *
   * \return
   *    Error code.
   */
  PPSTEE_ERROR create(int vertlocnbr,
                      int* vertloctab,
                      int* edgeloctab
                      );

  ///@}
  //! \name Getters/Setters
  ///@{

  ///@}

//private:

  //! \name Private data fields
  ///@{

  ///@}
};

#define PPSTEE_NUM_GID_ENTRIES 1  //!< Number of global index entries (for Zoltan).
#define PPSTEE_NUM_LID_ENTRIES 1  //!< Number of local index entries (for Zoltan).

//! PPStee Zoltan graph class
/*!
 * Graph class for access in native Zoltan description.
 */
class PPSteeGraphZoltan : public PPSteeGraph {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Full-parameter constructor. Calls create().
  /*!
   * Constructs a graph using the provided functions of Zoltan query function types.
   *
   * If #dataAccess is #PPSTEE_DATAACCESS_COPY functions are queried immediately and the obtained graph data is saved in the default arrays.
   * The query functions can be freed upon return.\n
   * If #dataAccess is #PPSTEE_DATAACCESS_QRYFN (references to) query functions are saved for later reference.\n
   * Setting #dataAccess \em directly to #PPSTEE_DATAACCESS_VIEW is \em not supported;
   * #dataAccess will be set to #PPSTEE_DATAACCESS_QRYFN.
   * (However, #dataAccess==#PPSTEE_DATAACCESS_VIEW is possible if the PPSteeGraphZoltan was derived from another graph type with #dataAccess==#PPSTEE_DATAACCESS_VIEW.
   *
   * \param[in] mpiComm MPI communicator.
   * \param[in] zNumObjFn, zObjListFn, zNumEdgesMultiFn, zEdgeListMultiFn zoltan query functions.
   * \param[in] dNumObj, dObjList, dNumEdgesMulti, dEdgeListMulti Pointers to user-defined data. May be omitted.
   * \param[in] access Type of data access.
   *
   */
  PPSteeGraphZoltan(MPI_Comm mpiComm,
      ZOLTAN_NUM_OBJ_FN * zNumObjFn,
      ZOLTAN_OBJ_LIST_FN* zObjListFn,
      ZOLTAN_NUM_EDGES_MULTI_FN* zNumEdgesMultiFn,
      ZOLTAN_EDGE_LIST_MULTI_FN* zEdgeListMultiFn,
      void* dNumObj = NULL,
      void* dObjList = NULL,
      void* dNumEdgesMulti = NULL,
      void* dEdgeListMulti = NULL,
      PPSTEE_DATAACCESS access = PPSTEE_DATAACCESS_QRYFN
      );

  //! Convert constructor.
  /*!
   * Construct Zoltan graph based on another graph. Especially, this includes graphs of other types.
   *
   * Conversion to Zoltan: uses Zoltan query functions. See http://www.cs.sandia.gov/Zoltan/ug_html/ug_query_lb.html
   *
   * \param[in] existingGraph
   */
  PPSteeGraphZoltan(PPSteeGraph* existingGraph);


  ///@}
  //! \name Graph operations.
  ///@{

  //! Create Zoltan graph.
  /*!
   *
   * \param[in] zNumObjFn, zObjListFn, zNumEdgesMultiFn, zEdgeListMultiFn zoltan query functions.
   * \param[in] dNumObj, dObjList, dNumEdgesMulti, dEdgeListMulti Pointers to user-defined data. May be omitted.
   *
   * \return
   *    Error code.
   */
  PPSTEE_ERROR create(
      ZOLTAN_NUM_OBJ_FN* zNumObjFn,
      ZOLTAN_OBJ_LIST_FN* zObjListFn,
      ZOLTAN_NUM_EDGES_MULTI_FN* zNumEdgesMultiFn,
      ZOLTAN_EDGE_LIST_MULTI_FN* zEdgeListMultiFn,
      void* dNumObj = NULL,
      void* dObjList = NULL,
      void* dNumEdgesMulti = NULL,
      void* dEdgeListMulti = NULL
      );

  ///@}
  //! \name Getters/Setters
  ///@{

  //! Get #vertglbtab.
  int* getVertglbtab();

  //! Get #vertloctab.
  int* getVertloctab();
  //! Get #edgeloctab.
  int* getEdgeloctab();


  ///@}

private:

  //! \name Data management functions
  ///@{
  

  /**
   * \brief Uses query functions to create copies of the graph data.
   *
   * @return Error code.
   */
  PPSTEE_ERROR createDataCopiesFromQueryFunctions();


  /**
   * \brief Queries provided Zoltan functions for local graph data.
   *
   * Uses #zNumObjFn for #vertloccnt, #zNumEdgesMultiFn for #edgeloccnt, #zNumEdgesMultiFn for #vertloctab and #zEdgeListMultiFn for #edgeloctab (cf. \ref doc_conversion_from_parmetis_like_to_zoltan).
   * The global part of the graph data, #vertglbtab, can not be obtained directly by calls to the query functions. (The first access to #vertglbtab assembles it automatically using #vertloccnt and MPI_Allgather; cf. PPSteeGraph::buildGlobalDataFromLocalData().)
   *
   * @return Error code.
   */
  PPSTEE_ERROR queryZoltanFunctionsForLocalGraphData();


  ///@}
  //! \name Private data fields
  ///@{
  

  ZOLTAN_NUM_OBJ_FN * zNumObjFn; //!< Reference to Zoltan query function NumObjFn.
  ZOLTAN_OBJ_LIST_FN * zObjListFn; //!< Reference to Zoltan query function ObjListFn.
  ZOLTAN_NUM_EDGES_MULTI_FN * zNumEdgesMultiFn; //!< Reference to Zoltan query function NumEdgesMultiFn.
  ZOLTAN_EDGE_LIST_MULTI_FN * zEdgeListMultiFn; //!< Reference to Zoltan query function EdgeListMultiFn.

  void* dNumObj; //!< Reference to user-defined data for Zoltan query function NumObj.
  void* dObjList; //!< Reference to user-defined data for Zoltan query function ObjList.
  void* dNumEdgesMulti; //!< Reference to user-defined data for Zoltan query function NumEdgesMulti.
  void* dEdgeListMulti; //!< Reference to user-defined data for Zoltan query function EdgeListMulti.


  ///@}
};


//! PPStee weights class
/*!
 * Weights are (/should be) bound to a graph and cover optional weights for vertices and edges of the graph.
 */
class PPSteeWeights {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Constructor.
  /*!
   * Default constructor. Weights should be based on a graph.
   *
   * \param[in] graph
   */
  PPSteeWeights(PPSteeGraph* graph);

  //! Destructor
  ~PPSteeWeights();


  ///@}
  //! \name Specify weights usage
  /*!
   * By default, vertex \em and edge weights are used. This can be reset to only one of these mutually exclusive choices and back.
   */
  ///@{


  //! Use only vertex weights.
  void useOnlyVertexWeights() {this->type = PPSTEE_WEIGHTS_ONLYVERTEX; };

  //! Use only edge weights.
  void useOnlyEdgeWeights() {this->type = PPSTEE_WEIGHTS_ONLYEDGE; };

  //! Use both, vertex \em and edge, weights.
  void useVertexAndEdgeWeights() {this->type = PPSTEE_WEIGHTS_ALL; };


  ///@}
  //! \name Weights data access
  ///@{

  //! Set weights data.
  /*!
   * \note Does \em not set #type. If necessary, use useVertexAndEdgeWeights().
   *
   * \param[in] vertexWeights,edgeWeights
   * \return Error code.
   */
  PPSTEE_ERROR setWeightsData(int* vertexWeights, int* edgeWeights);

  //! Set only vertex weights data.
  /*!
   * \note Does \em not set #type. If necessary, use useOnlyVertexWeights().
   *
   * \param[in] vertexWeights
   * \return Error code.
   */
  PPSTEE_ERROR setOnlyVertexWeightsData(int* vertexWeights);

  //! Set only edge weights data.
  /*!
   * \note Does \em not set #type. If necessary, use useEdgeWeights().
   *
   * \param[in] edgeWeights
   * \return Error code.
   */
  PPSTEE_ERROR setOnlyEdgeWeightsData(int* edgeWeights);


  ///@}
  //! \name Free data
  ///@{

  //! Frees memory of private arrays holding most of this class' data.
  /*!
   * \warning Use with care.
   */
  void freeMemoryIfHoldingDataCopies();


  ///@}
  //! \name Getters/Setters
  ///@{

  //! Get #type.
  PPSTEE_WEIGHTS getType() const {return type; };

  //! Get #dataAccess.
  PPSTEE_DATAACCESS getDataAccess() const {return dataAccess; };

  //! Get either weights' equipartitioned flag.
  bool areWeightsEquipartitioned(PPSTEE_GRAPHDATUM type) const { if (type == PPSTEE_GRAPHDATUM_VERTEX) { return vertexWeightsEquipartitioned; } else if (type == PPSTEE_GRAPHDATUM_EDGE) { return edgeWeightsEquipartitioned; } };
  //! Get #vertexWeightsEquipartitioned.
  bool areVertexWeightsEquipartitioned() const { return vertexWeightsEquipartitioned; };
  //! Get #edgeWeightsEquipartitioned.
  bool areEdgeWeightsEquipartitioned() const { return edgeWeightsEquipartitioned; };

  //! Get either weights.
  int* getWeights(PPSTEE_GRAPHDATUM type) const { if (type == PPSTEE_GRAPHDATUM_VERTEX) { return vertexWeights; } else if (type == PPSTEE_GRAPHDATUM_EDGE) { return edgeWeights; } };
  //! Get #vertexWeights.
  int* getVertexWeights() const {return vertexWeights; };
  //! Get #edgeWeights.
  int* getEdgeWeights() const {return edgeWeights; };

  //! Get #vertloccnt.
  int getVertloccnt() const {return vertloccnt; };
  //! Get #edgeloccnt.
  int getEdgeloccnt() const {return edgeloccnt; };

  ///@}

private:

  //! \name Private data fields
  ///@{

  PPSTEE_WEIGHTS type;  //!< Type of these weights.

  PPSTEE_DATAACCESS dataAccess; //!< Type of data access for these weights.

  int vertloccnt; //!< local vertex count. Inherited from constructing graph.
  int edgeloccnt;  //!< local edge count. Inherited from constructing graph.

  int* vertexWeights;  //!< Vertex weights.
  int* edgeWeights;  //!< Edge weights.

  bool vertexWeightsEquipartitioned;  //!< Specifies if vertex weights are equipartitioned, i.e., weights are uniformly distributed.
  bool edgeWeightsEquipartitioned;  //!< Specifies if edge weights are equipartitioned, i.e., weights are uniformly distributed.

  ///@}
  //! \name Compatibility check for weights to graph.
  ///@{

  PPSTEE_ERROR checkCompatibilityToGraph();


  ///@}


};

//! PPStee partition class
/*!
 * A partition is the result of a partitioning done by PPStee.
 *
 * \note Data access is limited to #PPSTEE_DATAACCESS_COPY.
 */
class PPSteePart {
public:

  //! \name Constructor/Destructor
  ///@{

  //! Constructor.
  PPSteePart(PPSteeGraph* graph, int* existingPartitioning);

  //! Constructor deriving partitioning from a graph.
  PPSteePart(PPSteeGraph* graph);

  //! Destructor
  ~PPSteePart();


  ///@}
  //! \name Partitioning data access
  ///@{

  //! Set partitioning data.
  /*!
   * \param[in] partData
   * \return Error code.
   */
  PPSTEE_ERROR setPartitioningData(int* partData);


  ///@}
  //! \name Free data
  ///@{

  //! Frees memory of private arrays holding most of this class' data.
  /*!
   * \warning Use with care.
   */
  void freeMemoryIfHoldingDataCopies();


  ///@}
  //! \name Getters/Setters
  ///@{

  //! Get #vertloccnt.
  int getVertloccnt() const {return vertloccnt; };

  //! Returns #partitioning array.
  int* getPartData() const {return partitioning; };

  ///@}

private:

  //! \name Private data fields
  ///@{

  int vertloccnt; //!< local vertex count. Inherited from constructing graph.

  int* partitioning;  //!< Partitioning array.

  ///@}
};


#endif /* PPSTEE_OBJECTS_HPP_ */
