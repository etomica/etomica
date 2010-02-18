package etomica.graph.test;

import java.util.HashMap;
import java.util.Map;

import static etomica.graph.model.Metadata.*;


import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.PartitionedIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.GraphIterator;
import etomica.graph.property.HasNoRootEdge;

public class PartitionedIteratorTest extends GraphIteratorTest {

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
    if (map == 1) {
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
  }

  public GraphIterator getColorIterator() {

    return new PartitionedIterator(coloredRootMap, coloredFieldMap);
  }

  public GraphIterator getMonoIterator() {

    if (MONO_VIA_PARTITIONS) {
      return new PartitionedIterator(monoRootMap, monoFieldMap);
    }
    return new DefaultIterator((byte) (monoRootMap.get(COLOR_CODE_0) + monoFieldMap.get(COLOR_CODE_0)),
        monoRootMap.get(COLOR_CODE_0));
  }

  /**
   * Partitions involving colored root and field nodes. These are the the new graphs where
   * both root and field nodes have distinguishing color amongst them. In other words, it
   * allows one to model any arbitrary distinguishing node properties as colors.
   */
  public void testCartesianColorGraphs1() {

    setupColorMap(1);
    // nodes = 2ra2fa1fb total of 3072 graphs after 698ms ==> no ground truth to check
    // against!
    // testNaive((byte) 5, getColorIterator());
  }

  public void testCartesianColorGraphs2() {

    setupColorMap(2);
    // nodes = 2ra1rb2fa1fb total of 294K graphs after 65s ==> no ground truth to check
    // against!
    // testNaive((byte) 6, getColorIterator());
  }

  public void testNoRootEdgeCartesianColorGraphs1() {

    setupColorMap(1);
    // nodes = 2ra2fa1fb total of 1536 graphs after 494ms ==> no ground truth to check
    // against!
    // testNaive((byte) 5, new PropertyFilter(getColorIterator(), new HasNoRootEdge()));
  }

  public void testNoRootEdgeCartesianColorGraphs2() {

    setupColorMap(2);
    // nodes = 2ra1rb2fa1fb total of 36864 graphs after 4s ==> no ground truth to check
    // against!
    // testNaive((byte) 6, new PropertyFilter(getColorIterator(), new HasNoRootEdge()));
  }

  public void testIsoFreeCartesianColorGraphs1() {

    setupColorMap(1);
    // nodes = 2ra2fa1fb total of 576 graphs after 1s ==> no ground truth to check
    // against!
    // testNaive((byte) 5, new IsomorphismFilter(getColorIterator()));
  }

  public void testIsoFreeCartesianColorGraphs2() {

    setupColorMap(2);
    // nodes = 2ra1rb2fa1fb total of ??? graphs after ???s (estimated 15min) ==> no ground
    // truth to check against!
    // testNaive((byte) 6, new IsomorphismFilter(getColorIterator()));
  }

  public void testIsoFreeNoRootEdgeCartesianColorGraphs1() {

    setupColorMap(1);
    // nodes = 2ra2fa1fb total of 288 graphs after 1s ==> no ground truth to check
    // against!
    // testNaive((byte) 5, new IsomorphismFilter(new PropertyFilter(getColorIterator(),
    // new HasNoRootEdge())));
  }

  public void testIsoFreeNoRootEdgeCartesianColorGraphs2() {

    setupColorMap(2);
    // nodes = 2ra1rb2fa1fb total of 6528 graphs after 2min ==> no ground truth to check
    // against!
    // testNaive((byte) 6, new IsomorphismFilter(new PropertyFilter(getColorIterator(),
    // new HasNoRootEdge())));
  }

  public void testNoRootEdgeIsoFreeCartesianColorGraphs1() {

    setupColorMap(1);
    // nodes = 2ra2fa1fb total of 288 graphs after 1s ==> no ground truth to check
    // against!
    // testNaive((byte) 5, new PropertyFilter(new IsomorphismFilter(getColorIterator()),
    // new HasNoRootEdge()));
  }

  public void testNoRootIsoFreeEdgeCartesianColorGraphs2() {

    setupColorMap(2);
    // nodes = 2ra1rb2fa1fb total of 6528 graphs after ???min (>> 2min) ==> no ground
    // truth to check against!
    // testNaive((byte) 6, new PropertyFilter(new IsomorphismFilter(getColorIterator()),
    // new HasNoRootEdge()));
  }

  /**
   * Partitions involving just one color. These are the classic graphs with root and field
   * nodes and no distinguishing color amongst them. In other words, it is the framework
   * supported by the virial Java applet. When generating mono graphs, use the default or
   * stored iterator. The partitioned approach is inherently more expensive and complex to
   * set up. In this test case, I use a flag to select which iterator to create for mono
   * graphs.
   */
  public void testCartesianMonoGraphs1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 296 graphs after 42ms
    // testNaive((byte) 5, getMonoIterator());
  }

  public void testCartesianMonoGraphs2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 32768 graphs after 2s (3s using partitions)
    // testNaive((byte) 6, getMonoIterator());
  }

  public void testNoRootEdgeCartesianMonoGraphs1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 512 graphs after 260ms
    // testNaive((byte) 5, new PropertyFilter(getMonoIterator(), new HasNoRootEdge()));
  }

  public void testNoRootEdgeCartesianMonoGraphs2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 4096 graphs after 1s
    // testNaive((byte) 6, new PropertyFilter(getMonoIterator(), new HasNoRootEdge()));
  }

  public void testIsoFreeCartesianMonoGraphs1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 240 graphs after 620ms
    // testNaive((byte) 5, new IsomorphismFilter(getMonoIterator()));
  }

  public void testIsoFreeCartesianMonoGraphs2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 6528 graphs after 65s
    // testNaive((byte) 6, new IsomorphismFilter(getMonoIterator()));
  }

  public void testIsoFreeNoRootEdgeCartesianMonoGraphs1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 120 graphs after 613ms
    // testNaive((byte) 5, new IsomorphismFilter(new PropertyFilter(getMonoIterator(), new
    // HasNoRootEdge())));
  }

  public void testIsoFreeNoRootEdgeCartesianMonoGraphs2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 816 graphs after 2s
    // testNaive((byte) 6, new IsomorphismFilter(new PropertyFilter(getMonoIterator(), new
    // HasNoRootEdge())));
  }

  public void testNoRootEdgeIsoFreeCartesianMonoGraphs1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 120 graphs after 648ms
    // testNaive((byte) 5, new PropertyFilter(new IsomorphismFilter(getMonoIterator()), new
    // HasNoRootEdge()));
  }

  public void testNoRootIsoFreeEdgeCartesianMonoGraphs2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 816 graphs after 81s ==> order of the filters matters!!!
    // testNaive((byte) 6, new PropertyFilter(new IsomorphismFilter(getMonoIterator()),
    // new HasNoRootEdge()));
  }
}