package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.operations.DropOrphanNodes;
import etomica.graph.property.IsDisconnected;

public class IsArticulationPoint {
	
	int[] list;
	IsDisconnected m = new IsDisconnected();

	public boolean isArticulationPoint(Graph g, int x){
		
		DropOrphanNodes dropper = new DropOrphanNodes();
		list = new int[10];
		Graph c = g.copy();
		//Deleting X
		for(int j = c.getOutDegree((byte) x)-1;	j>=0;	j--){
			//System.out.println("Deleting Edge of "+x+" with "+c.getOutNode((byte) x,(byte)  j));
			c.deleteEdge((byte) x,c.getOutNode((byte) x,(byte)  j));	//public void deleteEdge(byte fromNode, byte toNode);
		}
		
		c = dropper.apply(c);
		
		if(c.nodeCount()==0)return true;
		
	//	System.out.println("Articulation point dropped");
		
		for(int i = 0;(i<c.nodeCount()); i++){       // To find all bonds
			for(int j=0;j<c.getOutDegree( (byte) i);j++){
		//	System.out.println(i+","+c.getOutNode( (byte) i,(byte) j )+" are the new bonds after dropping orphan bonds");
			}
		}
		
		if(c.nodeCount()<g.nodeCount()-1) return true;
		else if(m.isDisconnected(c)) return true;
		else return false;
		
	/*	System.out.println("********************");
		
		for(int i = 0;(i<c.nodeCount()); i++){       // To find all bonds
			for(int j=0;j<c.getOutDegree( (byte) i);j++){
			System.out.println(i+","+c.getOutNode( (byte) i,(byte) j ));
			}
		}
		*/
	}	
}
