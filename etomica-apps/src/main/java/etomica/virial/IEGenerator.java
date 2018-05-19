/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;

public class IEGenerator {
    
	public static Graph relabel(Graph g){
	
		/*
         * swap node i and the first root node...
		 */
		
		if (g.getNode((byte)1).getType() == Metadata.TYPE_NODE_ROOT) return g.copy();

        byte otherrootnode=0;
		for(byte i=2;i<g.nodeCount();i++){
			if(g.getNode(i).getType()==Metadata.TYPE_NODE_ROOT) {
				otherrootnode = i;
				break;
			}
		}
		
		Relabel relabel = new Relabel();
		
		byte[] list = new byte[g.nodeCount()];
		for(byte i=0;i<g.nodeCount();i++){
		    list[i] = i;
		}
		list[1] = otherrootnode;
		list[otherrootnode] = 1;
				
		RelabelParameters params = new RelabelParameters(list);
		return relabel.apply(g, params);
	}
	
	public static Graph fourierMul(Graph a, Graph b, byte n) {
		
		MulFlexibleParameters mulFlexibleParameters = MulFlexibleParameters.makeParametersWithNodes(new char[0], (byte) 100, n, n);
		MulFlexible mulFlex = new MulFlexible();
		
		Graph product = mulFlex.apply(a, b, mulFlexibleParameters);

		product.getNode(n).setType(Metadata.TYPE_NODE_FIELD);
		return product;
	}

	/*
	 * Superimposes both nodes 0 and 1 from both graphs.
	 */
	public static Graph SuperImposeGraphs(Graph a,Graph b){
		
		byte newNodeCount = (byte)(a.nodeCount()+b.nodeCount()-2);
		
		Graph newGraph = GraphFactory.createGraph(newNodeCount);
		
		//copy graph 1
		for(byte i=0;i<a.nodeCount();i++){
			for(byte j=(byte)(i+1);j<a.nodeCount();j++){
				//System.out.println(i+ "  "+j+ "   "+newNodeCount);
				if(a.hasEdge(i, j)) {
				    newGraph.putEdge(i, j);
				    newGraph.getEdge(i,j).setColor(a.getEdge(i,j).getColor());
				}
			}
		}		
		
		byte[] oldlist = new byte[b.nodeCount()];
		byte[] newlist = new byte[b.nodeCount()];
		
		for(byte i=0;i<b.nodeCount();i++){
			oldlist[i]=i;
		}
		
		byte newid = a.nodeCount();
		newlist[0] = 0;
		newlist[1] = 1;
		for(byte i=2;i<b.nodeCount();i++){
			newlist[i] = newid;
			newid++;
		}
		
		for(int i=0;i<newlist.length;i++){
			for(int j=i+1;j<newlist.length;j++){
				if(b.hasEdge(oldlist[i], oldlist[j])) {
				    newGraph.putEdge(newlist[i], newlist[j]);
                    newGraph.getEdge(newlist[i], newlist[j]).setColor(b.getEdge(oldlist[i], oldlist[j]).getColor());
				}
			}
		}
		
		if( (a.hasEdge((byte)0, (byte)1)) || (b.hasEdge((byte)0, (byte)1)) ) {
		    newGraph.putEdge((byte)0, (byte)1);
		    char c = a.hasEdge((byte)0, (byte)1) ? a.getEdge((byte)0, (byte)1).getColor() : b.getEdge((byte)0, (byte)1).getColor();
		    newGraph.getEdge((byte)0, (byte)1).setColor(c);
		}
	
		newGraph.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		newGraph.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
	
		newGraph.coefficient().multiply(a.coefficient());
		newGraph.coefficient().multiply(b.coefficient());
		
		return newGraph;
	}

	public static Set<Graph> Multiply(Set<Graph> a,Set<Graph> b){ // Adds the bonds.. no "multiplication"
		
		Set<Graph> product = new GraphList();
				
		for(Graph g : a){
			for(Graph h : b){
		
				int newNodeCount=0;
				newNodeCount = g.nodeCount()+h.nodeCount()-2;
				if(newNodeCount<=8){// if its greater than 16 we get an error. NodeCOunt = 17 is not req unless we need B19.
					Graph newGraph = IEGenerator.SuperImposeGraphs(g.copy(), h.copy());
					if(newGraph.nodeCount()==0) throw new RuntimeException("The product of the graphs is empty");
					product.add(newGraph);
				}
				//else System.out.println("Too big");
			}
		}
		
		if(product.isEmpty()) throw new RuntimeException("The product of the graphs is empty");
		return product;
				
	}
	
}
