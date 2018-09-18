/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Set;

import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.DropOrphanNodes;
import etomica.graph.operations.GraphOp;
import etomica.graph.operations.MaxIsomorph;
import etomica.graph.operations.MaxIsomorph.MaxIsomorphParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.traversal.BCVisitor;
import etomica.graph.util.GraphNumber;
import etomica.graph.viewer.ClusterViewer;

/**
 * Class that pretends to be a property and determines how to decompose a
 * diagram into parts that can be treated with FFT and others that can not.
 *
 * @author Andrew Schultz
 */
public class FFTDecomposition implements Property {

  protected IsBiconnected isBi;
  protected List<Byte> articulationPair;
  protected final MaxIsomorph maxIsomorph;
  protected final MaxIsomorphParameters mip;
  protected final Relabel relabel;
  protected List strands;
  protected final List<Integer> segmentList;
  protected Graph graphSegment;
  protected Comparator comparator;

  public FFTDecomposition() {
    isBi = new IsBiconnected();
    maxIsomorph = new MaxIsomorph();
    mip = new MaxIsomorphParameters(new GraphOp.GraphOpNull(), new Property() {
      public boolean check(Graph graph) {
        return graph.getNode((byte)0).getType() == Metadata.TYPE_NODE_ROOT && graph.getNode((byte)1).getType() == Metadata.TYPE_NODE_ROOT;
      }
    });
    relabel = new Relabel();
    segmentList = new ArrayList<Integer>();
    comparator = new ComparatorStrand();
  }

  public List<Byte> getArticulationPair() {
    return articulationPair;
  }

  public List getStrands() {
    return strands;
  }

  public List<Integer> getSegments() {
    return segmentList;
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

  public int testPair(Graph g, byte start, byte end, List myStrands) {
    int maxSegment = 0;
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
    // make the copy here while we have start/end labeled as root
    graphSegment = g.copy();
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

    if (iStrand < 2) {
      if (oldPaths.size()+2 < g.nodeCount()) {
        // some of the nodes in g are not actually between start+end, delete them from graphSegment
        for (Node node : graphSegment.nodes()) {
          byte id = node.getId();
          if (id != start && id != end && !oldPaths.contains(id)) {
            for (byte i=0; i<graphSegment.getOutDegree(id); i++) {
              byte j = graphSegment.getOutNode(id, i);
              graphSegment.deleteEdge(id, j);
            }
          }
        }
        DropOrphanNodes don = new DropOrphanNodes();
        graphSegment = don.apply(graphSegment);
      }
      // graphSegment now has the bit we determined to be irreducible
      return Integer.MAX_VALUE;
    }
    
    HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
    
    for (byte k=0; k<visited.size(); k++) {
      // create a Graph for this strand
      List<Byte> kVisited = visited.get(k);
      if (kVisited.size() < 4) {
        myStrands.add(kVisited.size());
        if (kVisited.size() > maxSegment) maxSegment = kVisited.size();
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
      // gk should now have all the edges in the k strand

      hap.check(gk);
      List<Byte> aps = hap.getArticulationPoints();
      List<List<Byte>> biComponents = BCVisitor.getBiComponents(gk);
      List segments = new ArrayList();

      // we could factor gk, but we wouldn't know the articulation points that
      // connect the segments.  So just get the articulation points and use
      // them to define the segments
      for (List<Byte> lBiComponent : biComponents) {
        // each component should contain 2 points of interest,
        // articulation points or end points
        if (lBiComponent.size() < 3) {
          segments.add(lBiComponent.size());
          if (lBiComponent.size() > maxSegment) maxSegment = lBiComponent.size();
          continue;
        }
        byte lStart = -1, lEnd = -1;
        for (byte i : lBiComponent) {
          if (i == kStart || i == kEnd || aps.contains(i)) {
            if (lStart == -1) lStart = i;
            else lEnd = i;
          }
        }
        List lStrands = new ArrayList();
        int lMaxSegment = testPair(gk, lStart, lEnd, lStrands);
        if (lMaxSegment == Integer.MAX_VALUE) {
          // l was an irreducible segment, now stored as graphSegment
          graphSegment = maxIsomorph.apply(graphSegment, mip);
          int gn = Integer.parseInt(graphSegment.getStore().toNumberString());
          segments.add(gn);
          if (gn > maxSegment) maxSegment = gn;
          continue;
        }
        if (lMaxSegment > maxSegment) {
          // l contained a segment that was more complicated than our previous max
          maxSegment = lMaxSegment;
        }
        segments.add(lStrands);
      }
      if (segments.size() > 1) {
        Collections.sort(segments, comparator);
        myStrands.add(segments);
      }
      else {
        myStrands.add(segments.get(0));
      }
    }
    Collections.sort(myStrands, comparator);
    return maxSegment;
  }
  
  public boolean check(Graph graph) {

    if (graph.nodeCount() == 0) {
      return false;
    }
    if (!isBi.check(graph)) {
      return false;
    }
    if (graph.nodeCount() < 3) {
      return true;
    }
    
    byte nodeCount = graph.nodeCount();
    List ijStrands = new ArrayList();
    strands = null;
    int maxSegment = Integer.MAX_VALUE;
    int maxSegmentDepth = Integer.MAX_VALUE;
    byte iMin = -1, jMin = -1;
//    System.out.println("hi "+graph);
    for (byte i=0; i<nodeCount; i++) {
      for (byte j=(byte)(i+1); j<nodeCount; j++) {
        int ijMaxSegment = testPair(graph, i, j, ijStrands);
        if (ijMaxSegment < 4) {
          articulationPair = new ArrayList<Byte>();
          articulationPair.add(i);
          articulationPair.add(j);
          strands = ijStrands;
          return true;
        }
        if (ijMaxSegment < maxSegment) {
          iMin = i;
          jMin = j;
          maxSegment = ijMaxSegment;
          strands = ijStrands;
          ijStrands = new ArrayList();
          maxSegmentDepth = findSegmentDepth(strands, ijMaxSegment);
//          printStrands(strands, 1);
//          System.out.println("found1 "+maxSegment+" at depth "+maxSegmentDepth);
        }
        else if (ijMaxSegment == maxSegment) {
          int ijMaxSegmentDepth = findSegmentDepth(ijStrands, ijMaxSegment);
//          printStrands(ijStrands, 1);
//          System.out.println("found2 "+ijMaxSegment+" at depth "+ijMaxSegmentDepth);
          if (ijMaxSegmentDepth < maxSegmentDepth) {
            iMin = i;
            jMin = j;
            maxSegment = ijMaxSegment;
            strands = ijStrands;
            ijStrands = new ArrayList();
            maxSegmentDepth = ijMaxSegmentDepth;
          }
          else {
            ijStrands.clear();
          }
        }
        else {
          ijStrands.clear();
        }
      }
    }
    articulationPair = new ArrayList<Byte>();
    articulationPair.add(iMin);
    articulationPair.add(jMin);
    if (!segmentList.contains(maxSegment)) {
      segmentList.add(maxSegment);
    }
    return true;
  }
  
  public static int findSegmentDepth(List list, int segment) {
    int myMaxDepth = -1;
    for (Object obj : list) {
      if (obj instanceof Number) {
        if (myMaxDepth > -1) continue;
        if (((Integer)obj) == segment) myMaxDepth = 1;
      }
      else {
        int listMaxDepth = findSegmentDepth((List)obj, segment);
        if (listMaxDepth>-1 && listMaxDepth+1 > myMaxDepth) myMaxDepth = listMaxDepth+1;
      }
    }
    return myMaxDepth;
  }

  public static void printStrands(List list, int deep) {
    boolean lastNum = false;
    for (Object obj : list) {
      if (obj instanceof Number) {
        if (lastNum) System.out.print(",");
        System.out.print(obj);
        lastNum=true;
      }
      else {
        System.out.print("(");
        printStrands((List)obj, deep+1);
        System.out.print(")");
        lastNum = false;
      }
    }
    if (deep==1) System.out.println();
  }
  
  public static void main(String[] args) {
    int n = 5;
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
    
    FFTDecomposition isFFT = new FFTDecomposition();
    Set<Graph> fftSet = new GraphList(comp);
    System.out.println("FFT");
    for (Graph g : graphList) {
      if (!isFFT.check(g)) {
        // every diagram is FFT-able
        throw new RuntimeException("oops");
      }
      System.out.println(g);
      System.out.println("AP: "+isFFT.getArticulationPair());
      printStrands(isFFT.getStrands(), 1);
      fftSet.add(g);
    }
    ClusterViewer.createView("FFT", fftSet);
    Set<Graph> segmentGraphs = new GraphList(comp);
    for (long gn : isFFT.getSegments()) {
      Graph g = GraphNumber.makeGraph(gn);
      g.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      g.getNode((byte)1).setType(Metadata.TYPE_NODE_ROOT);
      segmentGraphs.add(g);
    }
    System.out.println("\nSegments:");
    for (Graph g : segmentGraphs) {
      System.out.println(g.getStore().toNumberString());
    }
    ClusterViewer.createView("segments", segmentGraphs);
  }
  
  /**
   * Comparator that recursively compares two strands
   */
  public static class ComparatorStrand implements Comparator {
    public int compare(Object o1, Object o2) {
      if (o1 instanceof Integer && o2 instanceof Integer) {
        int i1 = (Integer)o1;
        int i2 = (Integer)o2;
        return i1 > i2 ? 1 : (i1 == i2 ? 0 : -1); 
      }
      if (o1 instanceof Integer) {
        return -1;
      }
      else if (o2 instanceof Integer) {
        return 1;
      }
      List l1 = (List)o1;
      List l2 = (List)o2;
      int s1 = l1.size();
      int s2 = l2.size();
      for (int i=0; i<s1 && i<s2; i++) {
        int c = compare(l1.get(i), l2.get(i));
        if (c != 0) return c;
      }
      return s1 > s2 ? 1 : (s1 == s2 ? 0 : -1);
    }
  }
}