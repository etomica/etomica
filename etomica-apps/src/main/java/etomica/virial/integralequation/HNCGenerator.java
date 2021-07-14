/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.integralequation;

import etomica.graph.model.*;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.IsoFree;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

public class HNCGenerator {

	public static ArrayList<Set<Graph>> HNCGenerate(int n, char fBond) {
		
	
		Set<Graph> c0 = new GraphList();
		Set<Graph> t0 = new GraphList();
		Set<Graph> h0 = new GraphList();
		
//		MetadataImpl.rootPointLabelSpecial= false;
		MetadataImpl.rootPointsSpecial = false;
		        
		//create a Map
		ArrayList<Set<Graph>> cList = new ArrayList<Set<Graph>>();
		ArrayList<Set<Graph>> tList = new ArrayList<Set<Graph>>();
		ArrayList<Set<Graph>> hList = new ArrayList<Set<Graph>>();
		
		//fbond
		Graph fGraph = GraphFactory.createGraph((byte)2, BitmapFactory.createBitmap((byte)2,true));
		fGraph.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		fGraph.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
		fGraph.getEdge((byte)0, (byte)1).setColor(fBond);
				   
		//create a list  
		c0.add(fGraph);
		h0.add(fGraph);
		
		//populate the Map  
		cList.add(c0);
		tList.add(t0);
		hList.add(h0);
		
        Set<Graph> ttk = new GraphList();

        for(int m=1;m<=n;m++){

		    Set<Graph> cm = new HashSet<Graph>();
		    Set<Graph> tm = new HashSet<Graph>();
		    Set<Graph> hm = new HashSet<Graph>();

		//	System.out.println("m = "+m);
			
			//Calculating tn
			//R
		//	System.out.println("Calculating tn ");
			for(int i=0;i<m;i++){
				
				int j=m-i-1;
				
				Set<Graph> ci = cList.get(i);
				Set<Graph> hj = hList.get(j);
			//	ClusterViewer.createView("c"+i,cc);
				//ClusterViewer.createView("h"+j,hh);
			//	System.out.println("We are multiplying c"+i+" and h"+j+" to get t"+m);
				if(ci.isEmpty()) throw new RuntimeException(" cc is empty");
				if(hj.isEmpty()) throw new RuntimeException(" hh is empty");
				//int q=0;
				for(Graph gci:ci){
					for(Graph ghj:hj){
						//System.out.println("Graph #"+q);q++;
						
						if(gci.nodeCount()==0) throw new RuntimeException("gc is empty");
						if(ghj.nodeCount()==0) throw new RuntimeException("gh is empty");
						Graph newGraph = IEGenerator.fourierMul(gci.copy(),ghj.copy(),(byte)1);
					//	System.out.println("gc coeff = "+gc.coefficient()+" ;gh coeff = "+gh.coefficient()+"; newGraph = "+newGraph.coefficient());
						newGraph = IEGenerator.relabel(newGraph);

						tm.add(newGraph.copy());
						ttk.add(newGraph.copy()); // here k = 1
					}
				}
				
				//ClusterViewer.createView("t"+i+j,t);
					
			}
			
			tList.add(tm);
			
			//
			//*************************** Calculating hn************************************
			//
				
			
			if(m==1){
		//		System.out.println("Calculating h1");
				hm.addAll(ttk);
				Set<Graph> p = new GraphList(); // temp variable
				
				for(Graph g:hm){	
					Graph newGraph = g.copy();
					newGraph.putEdge((byte)0, (byte)1);// this is to take into account that hn = (f + 1)*(collection of ts)
					newGraph.getEdge((byte)0, (byte)1).setColor(fBond);
					p.add(newGraph);
				}
				
				for(Graph g:p){ // cant take graphs from h and add to h itself since it becomes an infinite loop.. so i used p
					hm.add(g);
				}
				
			}
			
			else {
	            Set<Graph> ttotal = new GraphList();
	            Set<Graph> tcount = new GraphList();

	            for(Graph g:ttk){  //Initialised to ttk
	                ttotal.add(g.copy());
	                tcount.add(g.copy());
	            }
					
				for(int i=2;i<=m;i++){
					
					tcount = IEGenerator.Multiply(tcount,ttk); 

					for(Graph g:tcount){
						int den = g.coefficient().getDenominator();
						g.coefficient().setDenominator(den*i); // 1/m!
						ttotal.add(g.copy());
					}
													
				}

			    for(Graph g:ttotal){
					if(g.nodeCount()==m+2){
						hm.add(g.copy());
					}
				}
					
	            Set<Graph> temp = new GraphList();
				for(Graph g:hm){
					temp.add(g);
					g = g.copy();
					g.putEdge((byte)0,(byte)1);// this is to take into account that hn = (f + 1)*(collection of ts)
					g.getEdge((byte)0,(byte)1).setColor(fBond);
					temp.add(g);
				}
				hm = temp;
			}
            hList.add(hm);
				
			// Calculating Cn
			
			Set<Graph> h = hList.get(m);
			
			for(Graph g: h){
				cm.add(g.copy());
			}
			for(Graph g:tm){
				Graph tempg = g.copy();
				tempg.coefficient().setNumerator(-tempg.coefficient().getNumerator());
				cm.add(tempg);
			}
		
			cList.add(cm);
			
		}
		
		return cList;
	}

	public static void main(String[] args){
		
		int n=3;
			
        ArrayList<Set<Graph>> cmap = HNCGenerator.HNCGenerate(n, 'f');
        Set<Graph> c = cmap.get(n);
        
        IsoFree Isofree = new IsoFree();
        
//        MetadataImpl.rootPointLabelSpecial= false;
        MetadataImpl.rootPointsSpecial = false;
            
        Set<Graph> temp = Isofree.apply(c, null);
       
        int Bn = n+2;
        
        for(Graph g:temp){
            
            Graph tempq = g.copy();
            g.coefficient().setDenominator(tempq.coefficient().getDenominator()*Bn);
            g.coefficient().setNumerator(tempq.coefficient().getNumerator()*(-1));
            System.out.println(g);
        }
	}
}
