package etomica.graph.iterators;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
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
      FileInputStream fis = new FileInputStream(new File(url.toURI()));
      ZipInputStream zis = new ZipInputStream(fis);
      ZipEntry ze = zis.getNextEntry();
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
    for (byte edgeId = 0; edgeId < nautyGraph.length(); edgeId++) {
      if (nautyGraph.charAt(edgeId) == '1') {
        g.putEdge(edgeId);
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