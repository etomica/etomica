/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.GraphOp;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.traversal.BCVisitor;
import etomica.graph.viewer.ClusterViewer;

/**
 * Property class that determines if a diagram can be decomposed by FFT.  If it
 * can be decomposed, it also provides the way it decomposed the graph.
 *
 * @author Andrew Schultz
 */
public class IsFFT implements Property {

  private IsBiconnected isBi;
  private List<Byte> articulationPair;
  protected final MaxIsomorph maxIsomorph;
  protected final MaxIsomorphParameters mip;
  protected final Relabel relabel;
  protected List<Integer> strands;

  public IsFFT() {
    isBi = new IsBiconnected();
    maxIsomorph = new MaxIsomorph();
    mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), new Property() {
      
      public boolean check(Graph graph) {
        return graph.getNode((byte)0).getType() == Metadata.TYPE_NODE_ROOT && graph.getNode((byte)1).getType() == Metadata.TYPE_NODE_ROOT;
      }
    });
    relabel = new Relabel();
  }

  public List<Byte> getArticulationPair() {
    return articulationPair;
  }

  public List getStrands() {
    return strands;
  }

  // attempts to walk from start to end taking all possible paths.
  // start was arrived from reverseStart; if we encounter reverseStart, we
  // ignore it (refuse to visit).  nodes we visit are stored in visitedNodes.
  protected void walk(Graph g, List<Byte> visitedNodes, byte start, byte end, byte reverseStart) {
    for (byte k=0; k<g.getOutDegree(start); k++) {
      byte currentNode = g.getOutNode(start, k);
      if (currentNode == end) {
        if (!visitedNodes.contains(end)) {
          visitedNodes.add(end);
        }
        // just because we reached end doesn't mean we're done.
        // we need to look for other paths
        continue;
      }
      if (currentNode == reverseStart || visitedNodes.contains(currentNode)) continue;
      visitedNodes.add(currentNode);
      walk(g, visitedNodes, currentNode, end, start);
    }
  }

  protected boolean testPair(Graph g, byte start, byte end, List<Integer> myStrands) {
    g = g.copy();
    if (start != 0 || end != 1) {
      byte[] relabels = new byte[g.nodeCount()];
      for (byte i=0; i<g.nodeCount(); i++) {
        relabels[i] = i;
      }
      relabels[0] = start;
      relabels[start] = 0;
      if (start == 1) {
        relabels[0] = end;
      }
      else {
        relabels[1] = end;
      }
      relabels[end] = 1;
      g = relabel.apply(g, new RelabelParameters(relabels));
    }
    start = 0;
    end = 1;
    boolean oldValue = MetadataImpl.rootPointsSpecial;
    MetadataImpl.rootPointsSpecial = true;
    g.getNode(start).setType(Metadata.TYPE_NODE_ROOT);
    g.getNode(end).setType(Metadata.TYPE_NODE_ROOT);
    g = maxIsomorph.apply(g, mip);
    g.getNode(start).setType(Metadata.TYPE_NODE_FIELD);
    g.getNode(end).setType(Metadata.TYPE_NODE_FIELD);
    MetadataImpl.rootPointsSpecial = oldValue;
    
    List<List<Byte>> visited = new ArrayList<List<Byte>>();
    List<Byte> oldPaths = new ArrayList<Byte>();
    int iStrand = 0;
    // consider i,j as an articulation pair.
    for (byte k=0; k<g.getOutDegree(start); k++) {
      // consider leaving from start in all possible directions
      byte kNode = g.getOutNode(start,k);
      if (kNode == end) {
        iStrand++;
        myStrands.add(2);
        continue;
      }
      // ignore going to a node we've already seen
      if (oldPaths.contains(kNode)) continue;
      // kVisited is the nodes we see along the k path
      List<Byte> kVisited = new ArrayList<Byte>();
      kVisited.add(start);
      kVisited.add(kNode);
      walk(g, kVisited, kNode, end, start);
      if (!kVisited.contains(end)) {
        // was unable to get from start to end via kNode
        continue;
      }
      oldPaths.addAll(kVisited);
      oldPaths.remove((Byte)start);
      oldPaths.remove((Byte)end);
      visited.add(kVisited);
      iStrand++;
    }

    if (iStrand < 2) return false;
    
    HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
    
    for (byte k=0; k<visited.size(); k++) {
      // create a Graph for this strand
      List<Byte> kVisited = visited.get(k);
      if (kVisited.size() < 4) {
        myStrands.add(kVisited.size());
        continue;
      }
      Graph gk = GraphFactory.createGraph((byte)kVisited.size());
      byte kStart = -1, kEnd = -1;
      for (byte i=0; i<kVisited.size(); i++) {
        byte ki = kVisited.get(i);
        if (ki == start) {
          kStart = i;
        }
        else if (ki == end) {
          kEnd = i;
        }
        for (byte j=(byte)(i+1); j<kVisited.size(); j++) {
          byte kj = kVisited.get(j);
          if (kj == start) {
            kStart = j;
          }
          else if (kj == end) {
            kEnd = j;
          }
          if ((ki == start || ki == end) && (kj == start || kj == end)) continue;
          if (g.hasEdge(ki,kVisited.get(j))) {
            gk.putEdge(i,j);
          }
        }
      }

      hap.check(gk);
      List<Byte> aps = hap.getArticulationPoints();
      List<List<Byte>> biComponents = BCVisitor.getBiComponents(gk);
      List<Integer> segments = new ArrayList<>();
      
      for (List<Byte> lBiComponent : biComponents) {
        // each component should contain 2 points of interest,
        // articulation points or end points
        if (lBiComponent.size() < 3) {
          segments.add(lBiComponent.size());
          continue;
        }
        byte lStart = -1, lEnd = -1;
        for (byte i : lBiComponent) {
          if (i == kStart || i == kEnd || aps.contains(i)) {
            if (lStart == -1) lStart = i;
            else lEnd = i;
          }
        }
        List<Integer> lStrands = new ArrayList<>();
        if (!testPair(gk, lStart, lEnd, lStrands)) {
          return false;
        }
        segments.addAll(lStrands);
      }
      myStrands.addAll(segments);
    }  
    return true;
  }
  
  public boolean check(Graph graph) {

    if (graph.nodeCount() == 0) {
      return false;
    }
    if (!isBi.check(graph)) {
      return false;
    }
    if (graph.nodeCount() < 4) {
      return true;
    }
    
    byte nodeCount = graph.nodeCount();
    strands = new ArrayList<>();
    for (byte i=0; i<nodeCount; i++) {
      for (byte j=(byte)(i+1); j<nodeCount; j++) {
        boolean success = testPair(graph, i, j, strands);
        if (success) {
          articulationPair = new ArrayList<Byte>();
          articulationPair.add(i);
          articulationPair.add(j);
          return true;
        }
        strands.clear();
      }
    }
    return false;
  }
  
  public static void printStrands(List list, int deep) {
    for (Object obj : list) {
      if (obj instanceof Integer) {
//        for (int i=0; i<deep; i++) System.out.print(" ");
        System.out.print(obj);
      }
      else {
        System.out.print("(");
        printStrands((List)obj, deep+1);
        System.out.print(")");
      }
    }
    if (deep==1) System.out.println();
  }
  
  public static void main(String[] args) {
    int n = 6;
    if (args.length>0) {
      n = Integer.parseInt(args[0]);
    }
    
    IsBiconnected isBi = new IsBiconnected();
    GraphIterator iter = new PropertyFilter(new StoredIterator((byte)n), isBi);
    ComparatorChain comp = new ComparatorChain();
    comp.addComparator(new ComparatorNumFieldNodes());
    comp.addComparator(new ComparatorBiConnected());
    comp.addComparator(new ComparatorNumEdges());
    comp.addComparator(new ComparatorNumNodes());
    GraphList graphList = new GraphList(comp);
    MaxIsomorph maxIso = new MaxIsomorph();
    MaxIsomorphParameters mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), MaxIsomorph.PROPERTY_ALL);
    while (iter.hasNext()) {
      graphList.add(maxIso.apply(iter.next(), mip));
    }
    
    IsFFT isFFT = new IsFFT();
    Set<Graph> fftSet = new GraphList(comp);
    System.out.println("FFT");
    for (Graph g : graphList) {
      if (isFFT.check(g)) {
        System.out.println(g);
        System.out.println("AP: "+isFFT.getArticulationPair());
        printStrands(isFFT.getStrands(), 1);
        fftSet.add(g);
      }
    }
    ClusterViewer.createView("FFT", fftSet);

    System.out.println("not FFT");
    iter = new PropertyFilter(new StoredIterator((byte)n), isBi);
    Set<Graph> notFftSet = new GraphList(comp);
    while (iter.hasNext()) {
      Graph g = maxIso.apply(iter.next(), mip);
      if (!isFFT.check(g)) {
        System.out.println(g);
        notFftSet.add(g);
      }
    }
    ClusterViewer.createView("not FFT", notFftSet);
  }
}