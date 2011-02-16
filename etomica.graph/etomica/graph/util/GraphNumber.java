package etomica.graph.util;

import etomica.graph.model.Bitmap;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;


public class GraphNumber {

  public static void main(String[] args) {
    String usage = "usage: GraphNumber num [nodeCount]";
    if (args.length == 0) {
      System.out.println(usage);
      System.exit(1);
    }
    long graphNum = Long.parseLong(args[0]);
    Graph g;
    if (args.length == 2) {
      byte nodeCount = Byte.parseByte(args[1]);
      g = makeGraph(graphNum, nodeCount);
    }
    else {
      g = makeGraph(graphNum);
    }
    System.out.println(g.nodeCount()+" "+g.edgesToString());
  }

  public static Graph makeGraph(long graphNum) {
    int numBits = 64-Long.numberOfLeadingZeros(graphNum);
    byte nodeCount = 1;
    boolean success = false;
    for (nodeCount = 1; nodeCount<12; nodeCount++) {
      if (nodeCount*(nodeCount-1)/2 >= numBits) {
        success = true;
        break;
      }
    }
    if (!success) {
      throw new RuntimeException("invalid number");
    }
    return makeGraph(graphNum, nodeCount);
  }
    
  public static Graph makeGraph(long graphNum, byte nodeCount) {
    int numBits = 64-Long.numberOfLeadingZeros(graphNum);
    int nLeadingZeros = nodeCount*(nodeCount-1)/2 - numBits;
    String bitmapString = "";
    for (int i=0; i<nLeadingZeros; i++) {
      bitmapString += "0";
    }
    bitmapString += Long.toBinaryString(graphNum);
    Bitmap bitmap = BitmapFactory.createBitmap(bitmapString);
    return GraphFactory.createGraph(nodeCount, bitmap);
  }
}
