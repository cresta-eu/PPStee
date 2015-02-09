  The pre-processing interface PPStee is designed to balance the load of the overall simulation. It specifically includes all simulation parts, thus extends load balance from simulation core to pre-processing and post-processing tasks, and visualisation. The well-known third-party partitioning libraries ParMETIS, PTScotch and Zoltan are used to calculate the data distribution that is required for the load balance.

Basic usage:

Here, we describe how to use PPStee based on an example implementation. We assume an existing code that initialises its data and then does a standard ParMETIS call,

  ParMETIS_V3_PartKway(
      vtxdist, xadj, adjncy,
      vwgt, adjwgt,
      wgtflag, numflag, ncon, nparts,
      tpwgts, ubvec, options, edgecut,
      part,
      comm);

to retrieve a partitioning named "part". Other partitioners can be used analogously.
We start by initialising a "PPSteeGraph" object with the graph data we have, i.e. "vtxdist" for the global vertex distribution and "xadj" and "adjncy" for the thread-local adjacency structure:

  // get graph (as ParMETIS type)
  PPSteeGraph graph =
      PPSteeGraphParmetis(MPI_COMM_WORLD, vtxdist, xadj, adjncy);

Next, we construct weights objects derived from the graph as these have to be compatible. We fill in weights for the computation and visualisation part. These weights denote the work load (vertex weights, "xadj") and communication time (edge weights, "adjncy") each simulation part needs.

  // construct and set weights for computation
  PPSteeWeights wgtCmp(&graph);
  wgtCmp.setWeightsData(vwgt_c, adjwgt_c);
  // construct and set weights for visualisation
  PPSteeWeights wgtVis(&graph);
  wgtVis.setWeightsData(vwgt_v, adjwgt_v);

Now, we establish an instance of the interface’s main object and submit our graph and weights data.

  // get interface
  PPStee ppstee;
  // submit graph
  ppstee.submitGraph(graph);
  // submit weights
  ppstee.submitNewStage(wgtCmp, PPSTEE_STAGE_COMPUTATION);
  ppstee.submitNewStage(wgtVis, PPSTEE_STAGE_VISUALISATION);

Finally, we trigger the calculation of the partitioning and get the desired partitioning.

  // calculate partitioning
  PPSteePart* part;
  ppstee.getPartitioning(&part);


