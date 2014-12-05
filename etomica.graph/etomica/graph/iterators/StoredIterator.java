/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.iterators;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.jar.JarEntry;
import java.util.jar.JarFile;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphIterator;
import etomica.graph.model.impl.GraphImpl;

public class StoredIterator implements GraphIterator {

  private static final String GRAPH_BASE = "/etomica/graph/model/impl/";
  private static final String GRAPH_FILE = "graphset-";
  private static final String ZIP_FILE = "graph6fmt.zip";
  private byte data[] = null;
  private int ptr = -1;
  private byte nodeCount;
  private int coefficient;

  public StoredIterator(byte nodeCount) {

    this.nodeCount = nodeCount;
    final int BUFFER = 1024;
    assert (nodeCount >= 2 && nodeCount <= 9);
    coefficient = 1;
    for (int i = 1; i <= nodeCount; i++) {
      coefficient *= i;
    }
    URL url = getClass().getResource(GRAPH_BASE + ZIP_FILE);
    try {
      ZipInputStream zis = null;
      if (url.getProtocol().equals("jar")) {
      	URL myURL = new URL(url.getPath());
      	// myURL.getPath() will be something like
      	// jar:file:/usr/users/bob/ex_vir_hs.jar!/etomica/graph/model/impl/graph6fmt.zip
      	String jarFileName = myURL.getPath().split("!")[0];
      	JarFile jarFile = new JarFile(jarFileName);
      	// strip off leading slash
      	JarEntry jarEntry = jarFile.getJarEntry(myURL.getPath().split("!")[1].substring(1));
      	InputStream jis = jarFile.getInputStream(jarEntry);
        zis = new ZipInputStream(jis);
      }
      else {
        FileInputStream fis = new FileInputStream(new File(url.toURI()));
        zis = new ZipInputStream(fis);
      }
      ZipEntry ze;
      String gf = GRAPH_FILE + nodeCount + 'a';
      do {
        ze = zis.getNextEntry();
      } while (!ze.getName().equalsIgnoreCase(gf));
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

  public boolean hasNext() {

    return (ptr != -1);
  }

  public Graph next() {

    if (hasNext()) {
      // read the size of the automorphism group
      int autoGroupSize = readGroupSize();
      // second line: encoding of the graph is graph6 (check nauty manual);
      Graph g = toNativeGraph(readGraph6());
      g.coefficient().setNumerator(coefficient / autoGroupSize);
      // g is the representative of its automorphism group
      return g;
    }
    return null;
  }

  public void remove() {

    // no-op
  }

  private Graph toNativeGraph(String nautyGraph) {

    Graph g = new GraphImpl(nodeCount);
    byte toNode = 1;
    for (byte nautyEdgeId = 0; nautyEdgeId < nautyGraph.length(); nautyEdgeId++) {
      
      if (nautyGraph.charAt(nautyEdgeId) == '1') {
        byte fromNode = 0;
        
        for (; toNode < 15; toNode++) {
          // 0, 1, 3, 6, 10, 15, ...
          byte sectionStart = (byte) (toNode * (toNode - 1) / 2);
          // edge - candidate is the offset of the fromNode in a toNode section 
          if (nautyEdgeId < sectionStart + toNode) {
            fromNode = (byte)(nautyEdgeId - sectionStart);
            break;
          }
        } 
      
        g.putEdge(fromNode, toNode);
      }
    }
    // int index = 0;
    // for (byte value = 1; value <= (nodeCount - 1) + (nodeCount - 2); value++) {
    // for (byte x = 0; x < nodeCount; x++) {
    // for (byte y = (byte) (x + 1); y < nodeCount; y++) {
    // if (x + y == value) {
    // if (nautyGraph.charAt(index) == '1') {
    // g.putEdge((byte) x, (byte) y);
    // g.putEdge((byte) y, (byte) x);
    // }
    // index++;
    // }
    // }
    // }
    // }
    return g;
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
    assert (data[ptr] == 63 + nodeCount);
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
    int bitSize = nodeCount * (nodeCount - 1) / 2;
    if (graph.length() < bitSize) {
      graph += lower6Bits((byte) 0);
    }
    ptr = ++ptr2;
    if (ptr >= data.length) {
      ptr = -1;
    }
    return graph.substring(0, bitSize);
  }
}