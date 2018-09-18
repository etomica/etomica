/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

 package etomica.virial;

import java.util.Set;

import etomica.graph.model.Graph;
import etomica.graph.model.GraphList;
import etomica.graph.property.FactorFFT;
import etomica.graph.property.IsArticulationPair;
import etomica.graph.property.Property;
import etomica.virial.cluster.VirialDiagrams;

public class IsFFT implements Property {
	
	IsArticulationPair q = new IsArticulationPair();
	FactorFFT r = new FactorFFT();
	
	public boolean check(Graph g){

		for (int i =0;i<g.nodeCount();i++){  // Generate a list of Articulation Pairs
			for (int j=i+1;j<g.nodeCount();j++){
				if(q.isArticulationPair(g,i,j)){
					if(r.factor(g,(byte)i,(byte)j)) {
					    //	System.out.println("FFT possible with this AP");
						return true;
					}
				}
			}
		}
		return false;
	}
		
	public static void main(String[] args) {
		
		IsFFT m = new IsFFT();
		final int n =5;
        boolean multibody = false;
        boolean flex = false;
        Set<Graph> FFT,MSMC;
       	VirialDiagrams virialDiagrams = new VirialDiagrams(n, multibody, flex, false);
        virialDiagrams.setDoReeHoover(false);
        virialDiagrams.setDoShortcut(true);
        virialDiagrams.makeVirialDiagrams();
        Set<Graph> graphs = virialDiagrams.getVirialGraphs();
        FFT = new GraphList();
        MSMC = new GraphList();
       // ALL = new GraphList<Graph>();
        for( Graph g:graphs){
        	if((g.nodeCount()==n)){
        		//ALL.add(g);
        		if(m.check(g)){
        			//System.out.println("Hi");
        			FFT.add(g);
        			//System.out.println("FFT");
        			//System.out.println(g.getStore().toNumberString()+" "+g);
        		}
        		else {
        			//System.out.println("MSMC");
        			//System.out.println(g.getStore().toNumberString()+" "+g);
        			MSMC.add(g);
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

