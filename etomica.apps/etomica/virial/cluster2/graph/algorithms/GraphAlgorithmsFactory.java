package etomica.virial.cluster2.graph.algorithms;

import etomica.virial.cluster2.graph.algorithms.impl.BFTraversal;
import etomica.virial.cluster2.graph.algorithms.impl.CheckArticulatonPair;
import etomica.virial.cluster2.graph.algorithms.impl.CheckBiconnected;
import etomica.virial.cluster2.graph.algorithms.impl.CheckConnected;
import etomica.virial.cluster2.graph.algorithms.impl.CheckIsomorphic;
import etomica.virial.cluster2.graph.algorithms.impl.CheckNodalPoint;
import etomica.virial.cluster2.graph.algorithms.impl.DFTraversal;

/**
 * This factory maintains one instance of each algorithm. Since the algorithm
 * implementations have no state, this guarantees that only one instance of such
 * algorithms are ever created during program execution.
 * 
 * @author Demian Lessa
 * 
 */
public class GraphAlgorithmsFactory {

  private static GraphProperty algoArticulationPair = new CheckArticulatonPair();
  private static GraphProperty algoArticulationPoint = new CheckBiconnected();
  private static GraphProperty algoBiconnected = new CheckBiconnected();
  private static GraphProperty algoConnected = new CheckConnected();
  private static GraphProperty algoNodalPoint = new CheckNodalPoint();
  private static GraphPairProperty algoIsomorphic = new CheckIsomorphic();
  private static GraphTraversal algoBFT = new BFTraversal();
  private static GraphTraversal algoDFT = new DFTraversal();

  public static GraphProperty getArticulationPairAlgo() {

    return algoArticulationPair;
  }

  public static GraphProperty getArticulationPointAlgo() {

    return algoArticulationPoint;
  }

  public static GraphProperty getBiconnectedAlgo() {

    return algoBiconnected;
  }

  public static GraphProperty getConnectedAlgo() {

    return algoConnected;
  }

  public static GraphProperty getNodalPointAlgo() {

    return algoNodalPoint;
  }

  public static GraphPairProperty getIsomorphismAlgo() {

    return algoIsomorphic;
  }

  public static GraphTraversal getBFTAlgo() {

    return algoBFT;
  }

  public static GraphTraversal getDFTAlgo() {

    return algoDFT;
  }
}