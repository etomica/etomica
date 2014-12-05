/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.test;

import static etomica.graph.model.Metadata.COLOR_CODE_0;
import static etomica.graph.model.Metadata.COLOR_CODE_1;

import java.util.HashMap;
import java.util.Map;


import etomica.graph.iterators.DefaultIterator;
import etomica.graph.iterators.PartitionedIterator;
import etomica.graph.iterators.filters.IsomorphismFilter;
import etomica.graph.iterators.filters.PropertyFilter;
import etomica.graph.model.GraphIterator;
import etomica.graph.property.HasNoRootEdge;
import etomica.graph.property.IsBiconnected;
import etomica.graph.property.IsConnected;

public class ConnectivityTest extends GraphIteratorTest {

  protected void testNaive(byte nodeCount, GraphIterator iterator) {

    testTemplate(nodeCount, iterator);
  }

  private Map<Character, Byte> coloredFieldMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> monoFieldMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> coloredRootMap = new HashMap<Character, Byte>();
  private Map<Character, Byte> monoRootMap = new HashMap<Character, Byte>();

  public void reset() {

    super.reset();
    printPermutations = false;
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
    else if (map == 3) {
      monoRootMap.put(COLOR_CODE_0, (byte) 2);
      monoFieldMap.put(COLOR_CODE_0, (byte) 4);
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

    return new PartitionedIterator(monoRootMap, monoFieldMap);
  }

  public void testAllConnected() {

    reset();
    // nodes = 6: total of 26704 graphs after 364ms (virial applet takes several minutes)
    byte rangeBegin = 1;
    byte rangeEnd = 6;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
      // testNaive(i, new PropertyFilter(new PropertyFilter(new DefaultIterator(i), new
      // HasNoRootEdge()), new IsConnected()));
    }
  }

  public void testAllBiconnected() {

    reset();
    // nodes = 6: total of 11368 graphs after 557ms (virial applet takes several minutes)
    byte rangeBegin = 1;
    byte rangeEnd = 6;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
//      testNaive(i, new PropertyFilter(new PropertyFilter(new DefaultIterator(i), new HasNoRootEdge()),
//          new IsBiconnected()));
    }
  }

  public void testIsoFreeConnected() {

    reset();
    // nodes = 6: total of 112 graphs after 3s
    byte rangeBegin = 1;
    byte rangeEnd = 6;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
//      testNaive(i, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(new DefaultIterator(i),
//          new HasNoRootEdge()), new IsConnected())));
    }
  }

  public void testIsoFreeBiconnected() {

    reset();
    // nodes = 6: total of 56 graphs after 1s
    byte rangeBegin = 1;
    byte rangeEnd = 6;
    for (byte i = rangeBegin; i <= rangeEnd; i++) {
//      testNaive(i, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(new DefaultIterator(i),
//          new HasNoRootEdge()), new IsBiconnected())));
    }
  }

  /**
   * The examples below are graphs with both field and root nodes.
   */
  public void testAllConnectedWithRootNodes1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 314 graphs after 114ms
//    testNaive((byte) 5, new PropertyFilter(new PropertyFilter(getMonoIterator(), new
//    HasNoRootEdge()), new IsConnected()));
  }

  public void testAllConnectedWithRootNodes2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 2414 graphs after 743ms
//    testNaive((byte) 6, new PropertyFilter(new PropertyFilter(getMonoIterator(), new
//    HasNoRootEdge()), new IsConnected()));
  }

  public void testAllConnectedWithRootNodes3() {

    setupMonoMap(3);
    // nodes = 2r4f: total of 12424 graphs after 803ms
//    testNaive((byte) 6, new PropertyFilter(new PropertyFilter(getMonoIterator(), new
//    HasNoRootEdge()), new IsConnected()));
  }

  public void testAllBiconnectedWithRootNodes1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 74 graphs after 201ms
//    testNaive((byte) 5, new PropertyFilter(new PropertyFilter(getMonoIterator(), new
//    HasNoRootEdge()), new IsBiconnected()));
  }

  public void testAllBiconnectedWithRootNodes2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 404 graphs after 748ms
//    testNaive((byte) 6, new PropertyFilter(new PropertyFilter(getMonoIterator(), new
//    HasNoRootEdge()), new IsBiconnected()));
  }

  public void testAllBiconnectedWithRootNodes3() {

    setupMonoMap(3);
    // nodes = 2r4f: total of 4336 graphs after 873ms
//    testNaive((byte) 6, new PropertyFilter(new PropertyFilter(getMonoIterator(), new
//    HasNoRootEdge()), new IsBiconnected()));
  }

  public void testIsoFreeConnectedWithRootNodes1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 67 graphs after 548ms
//  testNaive((byte) 5, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(getMonoIterator(),
//  new HasNoRootEdge()), new IsConnected())));
  }

  public void testIsoFreeConnectedWithRootNodes2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 449 graphs after 1s
//  testNaive((byte) 6, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(getMonoIterator(),
//  new HasNoRootEdge()), new IsConnected())));
  }

  public void testIsoFreeConnectedWithRootNodes3() {

    setupMonoMap(3);
    // nodes = 2r4f: total of 701 graphs after 5s
//  testNaive((byte) 6, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(getMonoIterator(),
//  new HasNoRootEdge()), new IsConnected())));
  }

  public void testIsoFreeBiconnectedWithRootNodes1() {

    setupMonoMap(1);
    // nodes = 2r3f: total of 18 graphs after 244ms
//  testNaive((byte) 5, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(getMonoIterator(),
//  new HasNoRootEdge()), new IsBiconnected())));
  }

  public void testIsoFreeBiconnectedWithRootNodes2() {

    setupMonoMap(2);
    // nodes = 3r3f: total of 80 graphs after 1s
//  testNaive((byte) 6, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(getMonoIterator(),
//  new HasNoRootEdge()), new IsBiconnected())));
  }

  public void testIsoFreeBiconnectedWithRootNodes3() {

    setupMonoMap(3);
    // nodes = 2r4f: total of 251 graphs after 2s
//  testNaive((byte) 6, new IsomorphismFilter(new PropertyFilter(new PropertyFilter(getMonoIterator(),
//  new HasNoRootEdge()), new IsBiconnected())));
  }
}