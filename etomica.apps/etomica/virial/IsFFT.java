/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

 package etomica.virial;

import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.graph.property.FactorFFT;
import etomica.graph.property.IsArticulationPair;
import etomica.graph.property.IsArticulationPoint;
import etomica.graph.property.Property;
import etomica.graph.util.GraphNumber;
//import etomica.graph.viewer.ClusterViewer;
import etomica.virial.cluster.VirialDiagrams;

public class IsFFT implements Property {
	
	IsArticulationPoint p = new IsArticulationPoint();
	IsArticulationPair q = new IsArticulationPair();
	FactorFFT r = new FactorFFT();
	
	public boolean check(Graph g){
		
		int flag=0;
		Graph copy = g.copy();
		for (int i =0;((i<g.nodeCount())&&(flag==0));i++){  // Generate a list of Articulation Pairs
			for (int j=i+1;((j<g.nodeCount())&&(flag==0));j++){
				if((q.isArticulationPair(g,i,j))&&(flag==0)){
					//System.out.println("HELLO");
					g = copy.copy();
					flag=1;
				}
			}
		}
		
		if(flag==0) return false;
		flag=0;
		for (int i =0;((i<g.nodeCount())&&(flag==0));i++){  // Generate a list of Articulation Pairs
			for (int j=i+1;((j<g.nodeCount())&&(flag==0));j++){
				if((q.isArticulationPair(g,i,j))&&(flag==0)){
					g = copy.copy();
					if(r.factor(copy,(byte)i,(byte)j)) {
					//	System.out.println("FFT possible with this AP");
						continue;
					}
					else flag=1;
				}
		
			}
		}
		if( flag==0) return true;
		else return false;
	}
		
	public static void main(String[] args) throws IOException{
		
		IsFFT m = new IsFFT();
		final int n =5;
		int Count=0,FFTCount=0;
        boolean multibody = false;
        boolean flex = false;
        Set<Graph> FFT,MSMC;
       	VirialDiagrams virialDiagrams = new VirialDiagrams(n, multibody, flex, false);
        virialDiagrams.setDoReeHoover(false);
        virialDiagrams.setDoShortcut(true);
        virialDiagrams.makeVirialDiagrams();
        Set<Graph> graphs = virialDiagrams.getVirialGraphs();
        FFT = new GraphList<Graph>();
        MSMC = new GraphList<Graph>();
       // ALL = new GraphList<Graph>();
        for( Graph g:graphs){
        	if((g.nodeCount()==n)){
        		//ALL.add(g);
        		if(m.check(g)){
        			//System.out.println("Hi");
        			FFT.add(g);
        			//System.out.println("FFT");
        			//System.out.println(g.getStore().toNumberString()+" "+g);
        			FFTCount++;
        		}
        		else {
        			//System.out.println("MSMC");
        			//System.out.println(g.getStore().toNumberString()+" "+g);
        			MSMC.add(g);
        			Count++;
        		}
        	}
        }
        //ClusterViewer.createView("FFT", FFT);
      //  ClusterViewer.createView("MSMC", MSMC);
		//System.out.println(Count+" = MSMCCount");
		//System.out.println(FFTCount+" = FFTCount");
		
		/*int p=2092113;
		Graph g = GraphNumber.makeGraph(p,(byte) 7);
		if(m.check(g))System.out.println(p+ "--->FFT possible");else System.out.println(p+ "--->Not Possible by FFT");*/
	}
}

