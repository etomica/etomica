/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.operations.DropOrphanNodes;

public class IsArticulationPair {
	
	IsConnected m = new IsConnected();
	
	public boolean isArticulationPair(Graph g,int x, int y){

		Graph c = g.copy();
		DropOrphanNodes dropper = new DropOrphanNodes();
		//Deleting X
	//	System.out.println("Deleting the bonds of "+x);
		for(int i = c.getOutDegree((byte) x)-1;	i>=0;	i--){
           c.deleteEdge((byte) x,c.getOutNode((byte) x,(byte)  i));	//public void deleteEdge(byte fromNode, byte toNode);
		}
	
		//Deleting Y
	//	System.out.println("Deleting the bonds of "+y);
		for(int j = c.getOutDegree((byte) y)-1;	j>=0;	j--){
			if(c.hasEdge((byte) y,c.getOutNode((byte) y,(byte)  j)))	//public boolean hasEdge(byte fromNode, byte toNode);
				c.deleteEdge((byte) y,c.getOutNode((byte) y,(byte)  j));	//public void deleteEdge(byte fromNode, byte toNode);
		}
	
		c = dropper.apply(c);
		
		if(c.nodeCount()==0)return true;
	
		if(c.nodeCount()<g.nodeCount()-2) {
	//			System.out.println("c.nodeCount()<g.nodeCount()-2");
			return true;
		}
		if(!m.check(c)) {
	//		System.out.println("isDisconnected(c) is true");
			return true;
		}
	//		System.out.println("isDisconnected(c) is not true");
		return false;
	}

}
