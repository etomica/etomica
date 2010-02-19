package etomica.graph.test;

import static etomica.graph.model.Metadata.COLOR_CODE_0;
import static etomica.graph.model.Metadata.COLOR_CODE_1;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.IsomorphismPrefilteredPartitionedIterator;
import etomica.graph.iterators.IteratorToSet;
import etomica.graph.iterators.PartitionedIterator;
import etomica.graph.iterators.StoredIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.operations.Binary;
import etomica.graph.operations.Delete;


public class DeleteTest extends GraphIteratorTest {

  private static boolean MONO_VIA_PARTITIONS = false;
  private Map<Character, Byte> coloredFieldMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> monoFieldMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> coloredRootMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> monoRootMap = new HashMap<Character, Byte>();

  protected void testNaive(byte nodeCount, GraphIterator iterator) {

    testTemplate(nodeCount, iterator);
  }

  public void reset() {

    super.reset();
    printPermutations = true;
    printMemory = true;
    checkAssertion = false;

    coloredRootMap.clear();
    coloredFieldMap.clear();

    monoRootMap.clear();
    monoFieldMap.clear();
  }

  protected void setupMonoMap(int map) {

    reset();
    if (map == 1) {
      monoRootMap.put(COLOR_CODE_0, (byte) 2);
      monoFieldMap.put(COLOR_CODE_0, (byte) 3);
    }
    else if (map == 2) {
      monoRootMap.put(COLOR_CODE_0, (byte) 3);
      monoFieldMap.put(COLOR_CODE_0, (byte) 3);
    }
  }

  protected void setupColorMap(int map) {

    reset();
    if (map == 0) {
      coloredRootMap.put(COLOR_CODE_0, (byte) 0);
      coloredFieldMap.put(COLOR_CODE_0, (byte) 2);
      coloredFieldMap.put(COLOR_CODE_1, (byte) 1);
    }
    else if (map == 1) {
      coloredRootMap.put(COLOR_CODE_0, (byte) 2);
      coloredFieldMap.put(COLOR_CODE_0, (byte) 2);
      coloredFieldMap.put(COLOR_CODE_1, (byte) 1);
    }
    else if (map == 2) {
      coloredRootMap.put(COLOR_CODE_0, (byte) 2);
      coloredRootMap.put(COLOR_CODE_1, (byte) 1);
      coloredFieldMap.put(COLOR_CODE_0, (byte) 2);
      coloredFieldMap.put(COLOR_CODE_1, (byte) 1);
    }
    else if (map == 3) {
      coloredRootMap.put(COLOR_CODE_0, (byte) 1);
      coloredFieldMap.put(COLOR_CODE_0, (byte) 2);
      coloredFieldMap.put(COLOR_CODE_1, (byte) 1);
    }
  }

  public void testDelete1() {

    byte rangeBegin = 2;
    byte rangeEnd = 6;
    Binary delete = new Delete();
    IteratorToSet its = new IteratorToSet();
//    for (byte nodeCount = rangeBegin; nodeCount <= rangeEnd; nodeCount++) {
//      GraphIterator i1 = new IsomorphismFilter(new DefaultIterator(nodeCount));
//      GraphIterator i2 = new StoredIterator(nodeCount);
//      Set<Graph> s1 = its.getSet(i1);
//      Set<Graph> s2 = its.getSet(i2);
//      Set<Graph> r = delete.apply(s1, s2, null);
//      assertTrue(r.size() == 0);
//      show(s1, s2, r);
//      r = delete.apply(s2, s1, null);
//      assertTrue(r.size() == 0);
//      show(s2, s1, r);
//      r = delete.apply(s1, s1, null);
//      assertTrue(r.size() == 0);
//      show(s1, s1, r);
//      r = delete.apply(s2, s2, null);
//      assertTrue(r.size() == 0);
//      show(s2, s2, r);
//    }
  }

  public void testDelete2() {

    byte rangeBegin = 2;
    byte rangeEnd = 6;
    Binary delete = new Delete();
    IteratorToSet its = new IteratorToSet();
//    for (byte nodeCount = rangeBegin; nodeCount <= rangeEnd; nodeCount++) {
//      GraphIterator i1 = new DefaultIterator(nodeCount);
//      GraphIterator i2 = new StoredIterator(nodeCount);
//      Set<Graph> s1 = its.getSet(i1);
//      Set<Graph> s2 = its.getSet(i2);
//      Set<Graph> r = delete.apply(s2, s1, null);
//      assertTrue(r.size() == 0);
//      show(s2, s1, r);
//      r = delete.apply(s1, s1, null);
//      assertTrue(r.size() == 0);
//      show(s1, s1, r);
//      r = delete.apply(s2, s2, null);
//      assertTrue(r.size() == 0);
//      show(s2, s2, r);
//    }
  }

  public void testDelete3() {

//    Binary delete = new Delete();
//    IteratorToSet its = new IteratorToSet();
//    setupColorMap(0);
//    GraphIterator i1 = new PartitionedIterator(coloredRootMap, coloredFieldMap);
//    GraphIterator i2 = new IsomorphismPrefilteredPartitionedIterator(coloredRootMap, coloredFieldMap);
//    GraphIterator i3 = new IsomorphismFilter(new PartitionedIterator(coloredRootMap, coloredFieldMap));
//    GraphIterator i4 = new IsomorphismFilter(new IsomorphismPrefilteredPartitionedIterator(coloredRootMap, coloredFieldMap));
//    Set<Graph> g1 = its.getSet(i1);
//    Set<Graph> g2 = its.getSet(i2);
//    Set<Graph> g3 = its.getSet(i3);
//    Set<Graph> g4 = its.getSet(i4);
//    Set<Graph> r;
//    r = delete.apply(g2, g1, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g3, g1, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g4, g1, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g3, g2, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g4, g2, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g4, g3, null);
//    assertTrue(r.size() == 0);
  }

  public void testDelete4() {

//    Binary delete = new Delete();
//    IteratorToSet its = new IteratorToSet();
//    setupColorMap(3);
//    GraphIterator i1 = new PartitionedIterator(coloredRootMap, coloredFieldMap);
//    GraphIterator i2 = new IsomorphismPrefilteredPartitionedIterator(coloredRootMap, coloredFieldMap);
//    GraphIterator i3 = new IsomorphismFilter(new PartitionedIterator(coloredRootMap, coloredFieldMap));
//    GraphIterator i4 = new IsomorphismFilter(new IsomorphismPrefilteredPartitionedIterator(coloredRootMap, coloredFieldMap));
//    Set<Graph> g1 = its.getSet(i1);
//    Set<Graph> g2 = its.getSet(i2);
//    Set<Graph> g3 = its.getSet(i3);
//    Set<Graph> g4 = its.getSet(i4);
//    Set<Graph> r;
//    r = delete.apply(g2, g1, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g3, g1, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g4, g1, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g3, g2, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g4, g2, null);
//    assertTrue(r.size() == 0);
//    r = delete.apply(g4, g3, null);
//    assertTrue(r.size() == 0);
  }

  private void show(Set<Graph> s1, Set<Graph> s2, Set<Graph> result) {

    for (Graph g : result) {
      System.out.println(g);
    }
    System.out.println();
    System.out.println("s1 total.....: " + s1.size() + " graphs");
    System.out.println("s2 total.....: " + s2.size() + " graphs");
    System.out.println("del total....: " + result.size() + " graphs\n");
  }
}