package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.operations.DropOrphanNodes;

public class IsArticulationPoint {
	
	IsConnected m = new IsConnected();

	public boolean isArticulationPoint(Graph g, int x){
		
		DropOrphanNodes dropper = new DropOrphanNodes();
		Graph c = g.copy();
		//Deleting X
		for(int j = c.getOutDegree((byte) x)-1;	j>=0;	j--){
			//System.out.println("Deleting Edge of "+x+" with "+c.getOutNode((byte) x,(byte)  j));
			c.deleteEdge((byte) x,c.getOutNode((byte) x,(byte)  j));	//public void deleteEdge(byte fromNode, byte toNode);
		}
		
		c = dropper.apply(c);
		
		if(c.nodeCount()==0)return true;

		if(c.nodeCount()<g.nodeCount()-1) return true;
		if(!m.check(c)) return true;
		return false;
	}	
}
