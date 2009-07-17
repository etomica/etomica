package etomica.virial.cluster2.graph.impl;

import java.io.*;
import java.net.*;
import java.util.zip.*;

import etomica.virial.cluster2.bitmap.Bitmap;
import etomica.virial.cluster2.graph.Edges;
import etomica.virial.cluster2.graph.EdgesFilter;
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
    final int BUFFER = 2048;
    assert (factory.getNodeCount() >= 2 && factory.getNodeCount() <= 9);
    factory = edgesFactory;
    URL url = getClass().getResource(ZIP_FILE);
    try {
      FileInputStream fis = new FileInputStream(new File(url.toURI()));
      ZipInputStream zis = new ZipInputStream(new BufferedInputStream(fis,
          BUFFER));
      ZipEntry ze = zis.getNextEntry();
      String gf = GRAPH_FILE + factory.getNodeCount();
      do {
        if (ze != null && ze.getName().equalsIgnoreCase(gf)) {
          data = new byte[(int) ze.getSize()];
          zis.read(data, 0, data.length);
        }
        else {
          ze = zis.getNextEntry();
        }
      }
      while (!ze.getName().equalsIgnoreCase(gf));
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
    // number of isomorphisms: N!/automorphism_group_size
    double coefficient = 1;
    int automorphismGroupSize = readGroupSize(data);
    if (automorphismGroupSize > 0) {
      for (int i = 1; i <= factory.getNodeCount(); i++) {
        coefficient *= i;
      }
      coefficient /= automorphismGroupSize;
    }
    // second line: encoding of the graph is graph6 (check nauty manual);
    // this graph is the representative of its automorphism group
    return GraphFactory.nautyEdges(factory.getRepresentation(readGraph6(data)),
        coefficient);
  }

  // check graph6s specs
  private Bitmap readGraph6(byte[] source) {

    assert (ptr != data.length && data[ptr] != '\n');
    // if 0 <= n <= 62, the first byte is 63+n
    assert (data[ptr++] == 63 + factory.getNodeCount());
    String line = "";
    int ptr2 = ptr;
    while (ptr2 != data.length && data[ptr2] != '\n') {
      ptr2++;
    }
    // the bytes in the range [ptr, ptr2] encode the upper 
    // triangle representation of the adjancency matrix for 
    // this graph
    for (int i = ptr; i< ptr2; i++) {
      line += data[ptr++];
    }
    return null;

  
//    #define BIAS6 63
//    #define MAXBYTE 126
//    #define SMALLN 62
//    #define TOPBIT6 32
//    #define C6MASK 63

//    #define ADDELEMENT0(setadd,pos)  ((setadd)[SETWD(pos)] |= BITT[SETBT(pos)])
//    #define DELELEMENT0(setadd,pos)  ((setadd)[SETWD(pos)] &= ~BITT[SETBT(pos)])
//    #define FLIPELEMENT0(setadd,pos) ((setadd)[SETWD(pos)] ^= BITT[SETBT(pos)])
//    #define ISELEMENT0(setadd,pos) (((setadd)[SETWD(pos)] & BITT[SETBT(pos)]) != 0)
//    #define EMPTYSET0(setadd,m) \
//        {setword *es; \
//        for (es = (setword*)(setadd)+(m); --es >= (setword*)(setadd);) *es=0;}
//    #define GRAPHROW0(g,v,m) ((set*)(g) + (long)(v)*(long)(m))

//  m = (n + WORDSIZE - 1) / WORDSIZE;
//  WORDSIZE = 32; // guess
//    if (s[0] != ':')       /* graph6 format */
//    {
//        k = 1;
//        for (j = 1; j < n; ++j)
//        {
//            gj = GRAPHROW(g,j,m);
//      
//            for (i = 0; i < j; ++i)
//            {
//                if (--k == 0)
//                {
//              k = 6;
//              x = *(p++) - BIAS6;
//                }
//        
//                if (x & TOPBIT6)
//                {
//              gi = GRAPHROW(g,i,m);
//              ADDELEMENT(gi,j);
//              ADDELEMENT(gj,i);
//                }
//                x <<= 1;
//            }
//        }
//    }
  }

  private int readGroupSize(byte[] source) {

    assert (ptr != data.length && data[ptr] != '\n');
    String line = "";
    while (ptr != data.length && data[ptr] != '\n') {
      line += data[ptr++];
    }
    return Integer.valueOf(line);
  }
}