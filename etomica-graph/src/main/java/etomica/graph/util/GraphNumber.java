/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.util;

import etomica.graph.model.Bitmap;
import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.viewer.ClusterViewer;


public class GraphNumber {

  public static void main(String[] args) {

      args = new String[]{"-display", "27412"};

    String usage = "usage: GraphNumber [-display] num [nodeCount]";
    if (args.length == 0) {
      System.out.println(usage);
      System.exit(1);
    }
    String graphNumStr = args[0];
    boolean display = false;
    if (args.length > 1 && args[0].equals("-display")) {
      display = true;
      graphNumStr = args[1];
    }
    byte nodeCount = -1;
    if ((!display && args.length == 2) || (display && args.length == 3)) {
      nodeCount = Byte.parseByte(display ? args[2] : args[1]);
    }
    GraphList list = display ? new GraphList(null) : null;
    if (graphNumStr.contains(",")) {
      String[] graphNumArray = graphNumStr.split(",");
      for (int i=0; i<graphNumArray.length; i++) {
        handleGraph(graphNumArray[i], nodeCount, list);
      }
    }
    else {
      handleGraph(graphNumStr, nodeCount, list);
    }
    if (display){
    	ClusterViewer.createView("Graph "+graphNumStr, list);
    }
  }
  
  public static void handleGraph(String graphStr, byte nodeCount, GraphList list) {

	  String graphNumStr=graphStr.replaceAll("[A-Za-z]", "");
	  String graphLetterStr=graphStr.replaceAll("[0-9]", "");

	  Integer graphNum = Integer.valueOf(graphNumStr);
	  Graph g = makeGraph(graphNum,nodeCount);
	  
	  if(graphLetterStr.length()>0){
		  if (g.nodeCount()!=graphLetterStr.length()){
			  throw new RuntimeException("improper value of graphLetter string!");
		  }
		  for (byte m = 0; m<g.nodeCount();m++){
			  g.getNode(m).setColor(graphLetterStr.charAt(m));
		  }
	  }

	  System.out.println(g.nodeCount()+" "+g.edgesToString());
	  if (list != null) {
		  list.add(g);
	  }
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
    if (nodeCount==-1) return makeGraph(graphNum);
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
