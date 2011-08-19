package etomica.virial.cluster;


import java.util.HashSet;
import java.util.Set;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.Parameters;
import etomica.graph.operations.Unary;

public class CombineABSite implements Unary {
	
	public Set<Graph> apply(Set<Graph> argument, Parameters params) {


	    assert (params==null);
	    Set<Graph> result = new HashSet<Graph>();
	    for (Graph g : argument) {
	    	Graph newGraph = apply(g);
	      result.add(newGraph);
	    }
	    return result;
	  }

	  public Graph apply(Graph g) {
	      char capFACBond = MetadataImpl.edgeColorPairs.get(0).get(0);//1st pair, 1st bond from 1st pair
	      char capFCABond = MetadataImpl.edgeColorPairs.get(0).get(1);//1st pair, 2nd bond from 1st pair
	      char capFBCBond = MetadataImpl.edgeColorPairs.get(1).get(0);//2nd pair, 1st bond from 2nd pair
	      char capFCBBond = MetadataImpl.edgeColorPairs.get(1).get(1);//2nd pair, 2nd bond from 2nd pair	
	      
	      char mCapFACBond = MetadataImpl.edgeColorPairs.get(2).get(0);//3rd pair, 1st bond from 3rd pair
	      char mCapFCABond = MetadataImpl.edgeColorPairs.get(2).get(1);//3rd pair, 2nd bond from 3rd pair
	      char mCapFBCBond = MetadataImpl.edgeColorPairs.get(3).get(0);//4th pair, 1st bond from 4th pair
	      char mCapFCBBond = MetadataImpl.edgeColorPairs.get(3).get(1);//4th pair, 2nd bond from 4th pair

  		g=g.copy();
          for (byte j = 0;j<g.nodeCount();j++){//looping over nodes
        		int bondedACount = 0;//from 0 to ..
        		int bondedBCount = 0;//from 0 to ..
        		boolean siteAFirst = false;
        		for (byte i=(byte)0; i <g.nodeCount();i++){
        			if(i==j||!g.hasEdge(j, i))continue;
        			Edge edge = g.getEdge(j, i); 
        			if ((edge.getColor()== capFACBond)||(edge.getColor()== mCapFACBond)){
        				if(bondedBCount == 0){
        					siteAFirst = true;
        				}
        				bondedACount++;
        			}
    				else if ((edge.getColor()== capFBCBond)||(edge.getColor()== mCapFBCBond)) {
    					bondedBCount++;
    				}
        		}
					if((bondedBCount>bondedACount)||(bondedACount>0&&bondedACount==bondedBCount&&!siteAFirst)){
						for(byte i=(byte)0; i <g.nodeCount();i++){
		        			if(i==j||!g.hasEdge(j, i))continue;
		        			Edge edge = g.getEdge(j, i);
		        			Edge reverseEdge = g.getEdge(i, j); 
		        			if ((edge.getColor()== capFBCBond)){
								edge.setColor(capFACBond);
								reverseEdge.setColor(capFCABond);
		        			}
		        			else if ((edge.getColor()== mCapFBCBond)){
								edge.setColor(mCapFACBond);
								reverseEdge.setColor(mCapFCABond);
		        			}
		        			else if ((edge.getColor()== capFACBond)){
								edge.setColor(capFBCBond);
								reverseEdge.setColor(capFCBBond);
		        			}
		        			else if ((edge.getColor()== mCapFACBond)){
								edge.setColor(mCapFBCBond);
								reverseEdge.setColor(mCapFCBBond);
		        			}
						}

					}

          }
		return g;
	}
}
