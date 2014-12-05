/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster2.graph.impl;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.bitmap.BitmapFactory;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
import etomica.virial.cluster2.graph.EdgesRepresentation;
import etomica.virial.cluster2.graph.EdgesRepresentationFactory;
import etomica.virial.cluster2.graph.GraphFactory;

public class StoredEdgesGenerator extends AbstractEdgesGenerator {

  private static final String TAG_STORED = "Stored";
  private static final String GRAPH_FILE = "graphset-";
  private static final String ZIP_FILE = "graph6fmt.zip";
  private EdgesRepresentationFactory factory;
  private byte data[] = null;
  private int ptr = -1;

  public StoredEdgesGenerator(EdgesRepresentationFactory edgesFactory,
      EdgesFilter filter) {

    super(filter);
    final int BUFFER = 1024;
    assert (factory.getNodeCount() >= 2 && factory.getNodeCount() <= 9);
    factory = edgesFactory;
    URL url = getClass().getResource(ZIP_FILE);
    try {
      FileInputStream fis = new FileInputStream(new File(url.toURI()));
      ZipInputStream zis = new ZipInputStream(fis);
      ZipEntry ze = zis.getNextEntry();
      String gf = GRAPH_FILE + factory.getNodeCount() + 'a';
      do {
        ze = zis.getNextEntry();
      }
      while (!ze.getName().equalsIgnoreCase(gf));
      if (ze != null && ze.getName().equalsIgnoreCase(gf)) {
        ptr = 0;
        data = new byte[(int) ze.getSize()];
        int pt = 0;
        int nread = 0;
        byte[] buf = new byte[BUFFER];
        while ((nread = zis.read(buf, 0, BUFFER)) > -1) {
          System.arraycopy(buf, 0, data, pt, nread);
          pt += nread;
        }
      }
      zis.closeEntry();
      zis.close();
      if (data == null) {
        throw new FileNotFoundException();
      }
    }
    catch (FileNotFoundException e) {
      e.printStackTrace();
      throw new RuntimeException("Graph file not found");
    }
    catch (URISyntaxException e) {
      e.printStackTrace();
      throw new RuntimeException("URI syntax invalid");
    }
    catch (IOException e) {
      e.printStackTrace();
      throw new RuntimeException("IO error while reading zip file");
    }
  }

  @Override
  protected String getTag() {

    return TAG_STORED;
  }

  @Override
  protected Edges push() {

    if (ptr == -1) {
      return null;
    }
    // coefficient is the size of the automorphism group
    int coefficient = 1;
    int autoGroupSize = readGroupSize();
    if (autoGroupSize > 0) {
      for (int i = 1; i <= factory.getNodeCount(); i++) {
        coefficient *= i;
      }
      coefficient /= autoGroupSize;
    }
    // second line: encoding of the graph is graph6 (check nauty manual);
    // this graph is the representative of its automorphism group
    String nautyGraph = readGraph6();
    return GraphFactory.nautyEdges(toNativeGraph(nautyGraph), GraphFactory.defaultCoefficient(coefficient));
  }

  private EdgesRepresentation toNativeGraph(String nautyGraph) {
    
    int nodeCount = factory.getNodeCount();
    Bitmap store = BitmapFactory.getBitmap(nodeCount * (nodeCount - 1) / 2, false);
    EdgesRepresentation rep = factory.getRepresentation(store);
    int index = 0;
    for (int value = 1; value <= (nodeCount - 1) + (nodeCount - 2); value++) {
      for (int x = 0; x < nodeCount; x++) {
        for (int y = x+1; y < nodeCount; y++) {
          if (x + y == value) {
            if (nautyGraph.charAt(index) == '1') {
              store.setBit(rep.getEdgeID(x, y));
              store.setBit(rep.getEdgeID(y, x));
            }
            index++;
          }
        }
      }
    }
    return rep;
  }
  
  private int readGroupSize() {

    assert (ptr != data.length && data[ptr] != '\n');
    String line = "";
    while (ptr != data.length && data[ptr] != '\n') {
      line += data[ptr++] - '0';
    }
    ptr++;
    return Integer.valueOf(line);
  }

  private String lower6Bits(byte byteValue) {

    String result = "";
    for (int i = 5; i >= 0; i--) {
      int test = (byteValue >> i) & 1;
      result = result + test;
    }
    return result;
  }

  // read one graph6 from the byte stream
  private String readGraph6() {

    assert (ptr < data.length);
    // verify that N(n) = 63 + n
    assert (data[ptr] == 63 + factory.getNodeCount());
    ptr++;
    // build R(x)
    int ptr2 = ptr;
    while (ptr2 != data.length && data[ptr2] != '\n') {
      ptr2++;
    }
    String graph = "";
    for (int i = ptr; i < ptr2; i++) {
      assert (data[i] >= 63 && data[i] <= 126);
      // every ASCII character encodes a 6-bit value
      graph = graph + lower6Bits((byte) (data[i] - 63));
    }
    int bitSize = factory.getNodeCount() * (factory.getNodeCount() - 1) / 2;
    if (graph.length() < bitSize) {
      graph += lower6Bits((byte)0);
    }
    ptr = ++ptr2;
    if (ptr >= data.length) {
      ptr = -1;
    }
    return graph.substring(0, bitSize);
  }
}