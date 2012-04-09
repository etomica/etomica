package etomica.virial;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Coefficient;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.impl.MetadataImpl;
import etomica.graph.operations.IsoFree;

public class HNCGenerator {

	public Set<Graph> NodeCountGraphs(int n,Set<Graph> t){
		
		Set<Graph> a = new GraphList<Graph>();
		
		if(t.isEmpty())System.out.println("nodecountgraph method: incoming graph is empty");
		
		for(Graph g:t){
			if(g.nodeCount()==n){
				a.add(g);
			}
		}
		
		if(a.isEmpty())System.out.println("nodecountgraph method: outgoing graph is empty");
		return a;
		
	}
	
	public Graph addCoefficients(Coefficient a, Coefficient b){  //works..checked
		
		int p = a.getNumerator();
		int q = a.getDenominator();
		
		int m = b.getNumerator();
		int n = b.getDenominator();
		
		int i = p*n + m*q;
		int j = q*n;
		
		Graph c = GraphFactory.createGraph((byte)3, BitmapFactory.createBitmap((byte)3,true));
		c.coefficient().setNumerator(i);
		c.coefficient().setDenominator(j);
		
		return c;
		
	}
	
	public int factorial(int n){//works..checked
		
		int product=1;
		
		if (n==0) return 1;
		
		for( int i=1;i<=n;i++){
			product = product*i;
		}
		return product;
		
	}
	// Swap nodes..
		public Graph SwapNodes(int m,int n,Graph a){
			
			int nodeCount = a.nodeCount();
			
			Graph b = GraphFactory.createGraph((byte)nodeCount, BitmapFactory.createBitmap((byte)nodeCount,false));
			
		/*	System.out.println("a.NodeCount = "+ a.nodeCount() );
			System.out.println("b.NodeCount = "+ b.nodeCount() );
			System.out.println("m = "+ m );
			System.out.println("n = "+ n );
			System.out.println(a);*/
			
			for(int i=0;i<a.getOutDegree((byte)m);i++){
							
				if(n!=a.getOutNode((byte)m,(byte)i)){
			//		System.out.println(n + "  " + a.getOutNode((byte)m,(byte)i));
					b.putEdge((byte)n, a.getOutNode((byte)m,(byte)i));
			//		System.out.println("placing bond between "+n + " , " + a.getOutNode((byte)m,(byte)i));
				}
				
			}
			
		//	System.out.println("Done with m");
			
			for(int i=0;i<a.getOutDegree((byte)n);i++){
							
				if(m!=a.getOutNode((byte)n,(byte)i)){
			//		System.out.println(m+"   "+a.getOutNode((byte)n,(byte)i));
					b.putEdge((byte)m, a.getOutNode((byte)n,(byte)i));
			//		System.out.println("placing bond between "+m + " , " + a.getOutNode((byte)n,(byte)i));
				}
				
			}
			
		//	System.out.println("Done with n");
			
			for(int i=0;i<a.nodeCount();i++){
				
				if((i!=m)&&(i!=n)){
					
					for(int j=0;j<a.getOutDegree((byte)i);j++){
						
						if((a.getOutNode((byte)i, (byte)j)!=m)&&(a.getOutNode((byte)i,(byte) j))!=n){
			//				System.out.println("i = " + i);
			//				System.out.println("a.getOutNode((byte)i, (byte)j) "+a.getOutNode((byte)i, (byte)j));
							if(i!=a.getOutNode((byte)i,(byte) j)){
			//					System.out.println("placing bond between "+i+" , "+a.getOutNode((byte)i,(byte) j));
								b.putEdge((byte)i,a.getOutNode((byte)i,(byte) j));
							}
						}
					}
				}
				
			}
			
//			System.out.println("Done with other nodes");
			
			if(a.hasEdge((byte)m, (byte)n)){
				
				b.putEdge((byte)m, (byte)n);
		//		System.out.println("placing bond between "+m+" , "+n);
				
			}
			
			//for rootnodes
			for(int i=0;i<a.nodeCount();i++){
				     if(i==m)b.getNode((byte)n).setType(a.getNode((byte)m).getType());
				else if(i==n)b.getNode((byte)m).setType(a.getNode((byte)n).getType());
				else b.getNode((byte)i).setType(a.getNode((byte)i).getType());
			}
					
			return b;
			
		}

	public ArrayList<Set<Graph>> HNCGenerate(int n){
		
	
		IEGenerator IEGenerator = new IEGenerator();
		HNCGenerator HNCGenerator = new HNCGenerator();
		
		Set<Graph> c = new GraphList<Graph>();
		Set<Graph> t = new GraphList<Graph>();
		Set<Graph> ttotal = new GraphList<Graph>();
		Set<Graph> ttk = new GraphList<Graph>();
		Set<Graph> term = new GraphList<Graph>();
		Set<Graph> tcount = new GraphList<Graph>();
		//Set<Graph> finalt = new GraphList<Graph>();
		Set<Graph> h = new GraphList<Graph>();
		
		MetadataImpl.rootPointLabelSpecial= false;
		MetadataImpl.rootPointsSpecial = false;
		        
		//Set<Graph> temp = Isofree.apply(c, null);
		    
		//create a Map
		ArrayList<Set<Graph>> cmap = new ArrayList<Set<Graph>>();
		ArrayList<Set<Graph>> tmap = new ArrayList<Set<Graph>>();
		ArrayList<Set<Graph>> hmap = new ArrayList<Set<Graph>>();
		ArrayList<Set<Graph>> chumma = new ArrayList<Set<Graph>>();
		
		/*byte[] list = {(byte) 0, (byte) 2, (byte) 1};
		RelabelParameters params = new RelabelParameters(list);*/
		
		//fbond
		Graph fbond = GraphFactory.createGraph((byte)2, BitmapFactory.createBitmap((byte)2,true));
		fbond.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		fbond.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
				   
		//create a list  
		c.add(fbond);
		t.add(null);
		h.add(fbond);
		
		//populate the Map  
		cmap.add(c);
		tmap.add(t);
		hmap.add(h);
		
		//Clearing the sets since it has been added
		c = new HashSet<Graph>();
		t = new HashSet<Graph>();
		h = new HashSet<Graph>();
		
		for(int m=1;m<=n;m++){
			
		//	System.out.println("m = "+m);
			
			//Calculating tn
			//R
		//	System.out.println("Calculating tn ");
			for(int i=0;i<m;i++){
				
				int j=m-i-1;
				
		//		System.out.println("i = "+i+"; j = "+j);
				
				if(j<0){
		//			System.out.println("j is negative"); 
					break;
				}
				
				Set<Graph> cc = cmap.get(i);
				Set<Graph> hh = hmap.get(j);
			//	ClusterViewer.createView("c"+i,cc);
				//ClusterViewer.createView("h"+j,hh);
			//	System.out.println("We are multiplying c"+i+" and h"+j+" to get t"+m);
				if(cc.isEmpty())System.out.println(" cc is empty");
				if(hh.isEmpty())System.out.println(" hh is empty");
				//int q=0;
				for(Graph gc:cc){
					for(Graph gh:hh){
						//System.out.println("Graph #"+q);q++;
						
						Set<Graph> check = new GraphList<Graph>();
						if(gc.nodeCount()==0)System.out.println("gc is empty");
						if(gh.nodeCount()==0)System.out.println("gh is empty");
						check.add(gc);
						check.add(gh);
						Graph newGraph = IEGenerator.fourierMul(gc.copy(),gh.copy(),1);
					//	System.out.println("gc coeff = "+gc.coefficient()+" ;gh coeff = "+gh.coefficient()+"; newGraph = "+newGraph.coefficient());
						newGraph = IEGenerator.Relabel(newGraph).copy();
						chumma.add(check);
				
						check = new HashSet<Graph>();
						
						t.add(newGraph.copy());
						ttotal.add(newGraph.copy());// has all ts...No need for rho in here..
						ttk.add(newGraph.copy()); // here k = 1
						term.add(newGraph.copy());
												
					}
				}
				
				//ClusterViewer.createView("t"+i+j,t);
					
			}
			
			tmap.add(t);
			//Clearing the sets since it has been added
			t = new HashSet<Graph>();
			
			//
			//*************************** Calculating hn************************************
			//
				
			Set<Graph> temp = new GraphList<Graph>();
			
			if(m==1){
		//		System.out.println("Calculating h1");
				h = HNCGenerator.NodeCountGraphs(m+2, ttk);
				Set<Graph> p = new GraphList<Graph>(); // temp variable
				
				for(Graph g:h){	
					Graph newGraph = g.copy();
					newGraph.putEdge((byte)0, (byte)1);// this is to take into account that hn = (f + 1)*(collection of ts)
					p.add(newGraph);
				}
				
				for(Graph g:p){ // cant take graphs from h and add to h itself since it becomes an infinite loop.. so i used p
					h.add(g);
				}
				
				hmap.add(h);
				
				h = new HashSet<Graph>();
			}
			
			else {
					temp = new HashSet<Graph>();
					//finalt = new HashSet();
					ttotal = new HashSet<Graph>();
					tcount = new HashSet<Graph>();
					
					for(Graph g:ttk){  //Initialised to ttk
						
						ttotal.add(g.copy());
						tcount.add(g.copy());
					}
				//	if(m==2)ClusterViewer.createView("ttk",ttk);
					
					System.out.println("m = " +m);
					for(int i=2;i<=m;i++){
						
						System.out.println("i"+i);
					    //if(m==3)ClusterViewer.createView("tcountB",tcount);
						tcount = IEGenerator.Multiply(tcount,ttk); 
						//if(m==3)ClusterViewer.createView("tcountA",tcount);
						for(Graph g:tcount){
							int den = g.coefficient().getDenominator();
							g.coefficient().setDenominator(den*i); // 1/m!
							ttotal.add(g.copy());
						}
														
					}

				    for(Graph g:ttotal){
						if(g.nodeCount()==m+2){
							h.add(g.copy());
						}
					}
					
					for(Graph g:h){
						temp.add(g.copy());
						g.putEdge((byte)0,(byte)1);// this is to take into account that hn = (f + 1)*(collection of ts)
						temp.add(g.copy());
					}
		
					
					hmap.add(temp);
					//Clearing the sets since it has been added
					temp = new HashSet<Graph>();
					h = new HashSet<Graph>();
					term = new HashSet<Graph>();
					tcount = new HashSet<Graph>();
					ttotal = new HashSet<Graph>();
					
			}
				
			// Calculating Cn
			
			h = hmap.get(m);
			t = tmap.get(m);
			
			for(Graph g: h){
				c.add(g.copy());
			}
			for(Graph g:t){
				Graph tempg = g.copy();
				tempg.coefficient().setNumerator(tempg.coefficient().getNumerator()*(-1));
				c.add(tempg);
			}
		
			cmap.add(c);
			
			c = new HashSet<Graph>();
			t = new HashSet<Graph>();
			h = new HashSet<Graph>();
			
		}
		
		return cmap;
	}
	
	public void tryout(int n){
		
		HNCGenerator HNCGenerator = new HNCGenerator();
		
		int p=n;
		
		System.out.println("n = " +n);
		ArrayList<Set<Graph>> cmap = HNCGenerator.HNCGenerate(n);
        Set<Graph> c = cmap.get(p);
		
        IsoFree Isofree = new IsoFree();
		
		MetadataImpl.rootPointLabelSpecial= false;
	    MetadataImpl.rootPointsSpecial = false;
	        
	   Set<Graph> temp = Isofree.apply(c, null);
	   
	    int Bn = n+2;
		
		int count =0;
	    for(Graph g:temp){
	    	
	    	Graph tempq = g.copy();
	    	g.coefficient().setDenominator(tempq.coefficient().getDenominator()*Bn);
	    	g.coefficient().setNumerator(tempq.coefficient().getNumerator()*(-1));
	    	System.out.println(g);
	       	count++;
	    }
	        
	    //ClusterViewer.createView("c"+p,temp);
	  //  ClusterViewer.createView("1",temp);
	    System.out.println("Count = "+count);
		

	}
	
	public static void main(String[] args){
		
		int n=3;
		
		HNCGenerator HNCGenerator = new HNCGenerator();
			
	    HNCGenerator.tryout(n);
		//HNCGenerator.HNCGenerate(n);
		
	}
}
