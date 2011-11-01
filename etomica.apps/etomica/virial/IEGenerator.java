package etomica.virial;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.operations.MulFlexible;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.operations.Relabel;
import etomica.graph.operations.RelabelParameters;
import etomica.graph.viewer.ClusterViewer;

public class IEGenerator {
	
	public void GraphDetails(Graph g){
		
		
		//System.out.println(g);
		
	}
	
	public Graph Relabel(Graph g){
	
		/*
		 * This method relabels nodes 1 and 2.
		 */
		
		int otherrootnode=0;
		for(int i=1;i<g.nodeCount();i++){
			if(g.getNode((byte)i).getType()==Metadata.TYPE_NODE_ROOT) {
				otherrootnode = i;
				i=100;// to leave the loop.. dont feel like using break.
			}
		}
		
		Relabel Relabel = new Relabel();
		
		byte[] list = new byte[g.nodeCount()];
		for(int i=0;i<g.nodeCount();i++){
			if(i==1)list[i] = (byte)otherrootnode;
			else if(i==otherrootnode)list[i] = (byte)1;
			else list[i] = (byte)i;
		}
				
		RelabelParameters params = new RelabelParameters(list);
		g=Relabel.apply(g, params).copy();
		
		return g;
		
	}
	
	public Graph fourierMul(Graph p,Graph q,int n) {
		
		Graph a = p.copy();
		
		Graph b = q.copy();
		MulFlexibleParameters MulFlexibleParameters = new MulFlexibleParameters(new char[0], (byte) 100,(byte) n,(byte) n);
		MulFlexible MulFlexible = new MulFlexible();
		Relabel Relabel = new Relabel();
		IEGenerator IEGenerator = new IEGenerator();
		
		byte nodes[] = new byte[4];
		
		//RelabelParameters params = new RelabelParameters(0,1,2);
		
		// Setting color to the 1 node in both the graph
		
		char color = 'o';
		a.getNode((byte) 1).setColor(color);
		b.getNode((byte) 1).setColor(color);
		
		Graph product = MulFlexible.apply(a, b, MulFlexibleParameters);
		
		for(int i=0;i<product.nodeCount();i++){
			if(color==product.getNode((byte)i).getColor()){ 
				product.getNode((byte) i).setType(Metadata.COLOR_CODE_0);
				product.getNode((byte) i).setType(Metadata.TYPE_NODE_FIELD);
			}
		}
		//Relabel.apply(product,params);
			
		//product.coefficient().multiply(p.coefficient());
		//product.coefficient().multiply(q.coefficient());
		
		for(int i=0;i<product.nodeCount();i++){
			if((product.getNode( (byte) i).getType()!= Metadata.TYPE_NODE_ROOT)&&((product.getNode( (byte) i).getType()!= Metadata.TYPE_NODE_FIELD))){
				product.getNode((byte) i).setType(Metadata.COLOR_CODE_0);
				product.getNode((byte) i).setType(Metadata.TYPE_NODE_FIELD);
				//System.out.print(i+" is not root or field");
				
			}
		}
		
		Graph newGraph = GraphFactory.createGraph((byte)product.nodeCount(), BitmapFactory.createBitmap((byte)product.nodeCount(),false));
		
		for(int i=0;i<product.nodeCount();i++){
			for(int j=0;j<product.nodeCount();j++){
				if(product.hasEdge((byte)i,(byte)j))newGraph.putEdge((byte)i,(byte)j);
			}
		}
		
		for(int i=0;i<product.nodeCount();i++){
			if(product.getNode( (byte) i).getType()== Metadata.TYPE_NODE_ROOT){
				newGraph.getNode((byte) i).setType(Metadata.TYPE_NODE_ROOT);
			}
			else newGraph.getNode((byte) i).setType(Metadata.TYPE_NODE_FIELD);
		}
		
		newGraph.coefficient().setNumerator(product.coefficient().getNumerator());
		newGraph.coefficient().setDenominator(product.coefficient().getDenominator());
		
		return newGraph;
		
	}
	
	
	public Graph SuperImposeGraphs(Graph a,Graph b){
		
		IEGenerator IEGenerator = new IEGenerator();
		int newNodeCount=0,node1=0,node2=0;
		newNodeCount = a.nodeCount()+b.nodeCount()-2;
		
		Graph newGraph = GraphFactory.createGraph((byte) newNodeCount, BitmapFactory.createBitmap((byte)newNodeCount,false));
		
		//copy graph 1
		for(int i=0;i<a.nodeCount();i++){
			for(int j=0;j<a.nodeCount();j++)
			{
				//System.out.println(i+ "  "+j+ "   "+newNodeCount);
				if(a.hasEdge((byte)i, (byte)j)) newGraph.putEdge((byte)i, (byte)j);
			}
		}		
		
		int[] oldlist = new int[b.nodeCount()];
		int[] newlist = new int[b.nodeCount()];
		
		for(int i=0;i<b.nodeCount();i++){
			oldlist[i]=i;
		}
		
		int newid = a.nodeCount();
		for(int i=0;i<b.nodeCount();i++){
			
			if( (i!=1)&&(i!=0)){
				newlist[i] = newid;
				newid++;
			}
			else newlist[i] = i;
			
		}
		
		for(int i=0;i<newlist.length;i++){
			for(int j=0;j<newlist.length;j++){
				if(b.hasEdge((byte)oldlist[i], (byte)oldlist[j])) newGraph.putEdge((byte)newlist[i], (byte)newlist[j]);
			}
		}
		
		if( (a.hasEdge((byte)0, (byte)1)) || (b.hasEdge((byte)0, (byte)1)) ) newGraph.putEdge((byte)0, (byte)1);
	
		newGraph.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		newGraph.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
	
		newGraph.coefficient().multiply(a.coefficient());
		newGraph.coefficient().multiply(b.coefficient());
		
		return newGraph;
	}
	public Graph realMul(Graph a,Graph b){
		
		int newNodeCount=0,node1=0,node2=0;
		newNodeCount = a.nodeCount()+b.nodeCount()-2;
		Graph newGraph = GraphFactory.createGraph((byte) newNodeCount, BitmapFactory.createBitmap((byte)newNodeCount,false));
		
		//Getting the root nodes
		int flag=0;
		for(int i=0;i<a.nodeCount();i++){
			if(a.getNode((byte)i).getType()==Metadata.TYPE_NODE_ROOT && (flag ==0)) { node1 = i; flag =1;}
			if(a.getNode((byte)i).getType()==Metadata.TYPE_NODE_ROOT && (flag !=0)) node2 = i;
		}
		for(int i=0;i<a.nodeCount();i++){
			for(int j=0;j<a.nodeCount();j++)
			{
				if(a.hasEdge((byte)i, (byte)j)) newGraph.putEdge((byte)i, (byte)j);
			}
		}		
		
		int[] oldlist = new int[b.nodeCount()];
		int[] newlist = new int[b.nodeCount()];
		int list = a.nodeCount();
		for(int i=0;i<b.nodeCount();i++){
			oldlist[i]=i;
		}
		
		for(int i=0;i<oldlist.length;i++){
			if((oldlist[i]!=node1)&&(oldlist[i]!=node2)) {
				newlist[i] = list;
				list++;
			}
			else newlist[i] = oldlist[i];
		}
		for(int i=0;i<oldlist.length;i++){
			for(int j=0;j<oldlist.length;j++){
				if(a.hasEdge((byte)oldlist[i], (byte)oldlist[j])) newGraph.putEdge((byte)newlist[i], (byte)newlist[j]);
			}
		}
		newGraph.getNode((byte) node1).setType(Metadata.TYPE_NODE_ROOT);
		newGraph.getNode((byte) node2).setType(Metadata.TYPE_NODE_ROOT);
		
		newGraph.coefficient().multiply(a.coefficient());
		newGraph.coefficient().multiply(b.coefficient());
		
		return newGraph;
		
	}
	public Set<Graph> Multiply(Set<Graph> a,Set<Graph> b){ // Adds the bonds.. no "multiplication"
		
		IEGenerator IEGenerator = new IEGenerator();
		Set<Graph> product = new GraphList<Graph>();
				
		for(Graph g : a){
			for(Graph h : b){
		
				int newNodeCount=0,node1=0,node2=0;
				newNodeCount = g.nodeCount()+h.nodeCount()-2;
				if(newNodeCount<=8){// if its greater than 16 we get an error. NodeCOunt = 17 is not req unless we need B19.
					Graph newGraph = IEGenerator.SuperImposeGraphs(g.copy(), h.copy());
					if(newGraph.nodeCount()==0)System.out.println("The product of the graphs is empty");
					product.add(newGraph.copy());
				}
				//else System.out.println("Too big");
			}
		}
		
		if(product.isEmpty())System.out.println("The product of the graphs is empty");
		return product;
				
	}
	public Set<Graph> SetDetails(Set<Graph> sg){
		
		IEGenerator IEGenerator = new IEGenerator();
		int count=0;
		for(Graph g:sg){
		///	System.out.println("# "+count);
			IEGenerator.GraphDetails(g);
			count++;
		}
	//	System.out.println(" END OF LIST ");
		return null;
		
	}
	
	public ArrayList PYGenerate(int n){
		
		IEGenerator IEGenerator = new IEGenerator();
		Relabel Relabel = new Relabel();
		
		
		
		Set<Graph> c = new GraphList<Graph>();
		Set<Graph> t = new GraphList<Graph>();
		Set<Graph> h = new GraphList<Graph>();
		//create a Map
		ArrayList cmap = new ArrayList();
		ArrayList tmap = new ArrayList();
		ArrayList hmap = new ArrayList();
		
		byte[] list = {(byte) 0, (byte) 2, (byte) 1};
		RelabelParameters params = new RelabelParameters(list);
		
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
		c = new HashSet();
		t = new HashSet();
		h = new HashSet();
		
		for(int m=1;m<=n;m++){ //Go from C0 to Cn
			
		//	System.out.println(" m = "+m);
			//Calculating tn
			//R
		//	System.out.println(" R ");
			for(int i=0;i<m;i++){
				
				int j=m-i-1;
				
				if(j<0)break;
				
				Set<Graph> cc = (Set<Graph>) cmap.get(i);
				Set<Graph> hh = (Set<Graph>) hmap.get(j);
				//System.out.println("We are multiplying c"+i+" and h"+j+" to get t"+m);
				if(cc.isEmpty())System.out.println(i+"cc is empty");
				if(hh.isEmpty())System.out.println(j+"hh is empty");
				int q=0;
				for(Graph gc:cc){
					for(Graph gh:hh){
						//System.out.println("Graph #"+q);q++;
						if(gc.nodeCount()==0)System.out.println("gc is empty");
						if(gh.nodeCount()==0)System.out.println("gh is empty");
						Graph newGraph = IEGenerator.fourierMul(gc,gh,1);
						newGraph = IEGenerator.Relabel(newGraph).copy();
						IEGenerator.GraphDetails(newGraph);
						t.add(newGraph);
						h.add(newGraph);
							
					}
				}
					
			}
			if(t.isEmpty())System.out.println("t set is empty");
			if(h.isEmpty())System.out.println("h set is empty");
			
			//Adding the Graph Set to tmap
			tmap.add(t);
			//Clearing the sets since it has been added
			t = new HashSet();
			
			//Calculating c
			//P
			//System.out.println(" P ");
			for(int i=0;i<m;i++){
			
				int j=m-i-1;		
				
				if(j<0)break;
				
				Set<Graph> cc = (Set<Graph>) cmap.get(i);
				Set<Graph> kk = (Set<Graph>) cmap.get(j); //kk also contains c
				//System.out.println("We are multiplying c"+i+" and c"+j+" to get c"+m);
				if(cc.isEmpty())System.out.println(i+" P - cc is empty");
				if(kk.isEmpty())System.out.println(i+" P - kk is empty");
				int q=0;
				for(Graph gc:cc){
					for(Graph gk:kk){
						
						//System.out.println("Graph #"+q);q++;
						Graph temp = IEGenerator.fourierMul(gc,gk,1);
						temp = IEGenerator.Relabel(temp).copy();
						//IEGenerator.GraphDetails(temp);
						//System.out.println("With the fbond");
						temp.putEdge((byte)0,(byte)1);
						//temp = IEGenerator.realMul(temp,fbond).copy();
						IEGenerator.GraphDetails(temp);
						c.add(temp);
						h.add(temp);
						
					}
				}
				
				
			}
			//Q
			//System.out.println(" Q ");
			for(int i=0;i<m;i++){
			
				int j=m-i-1;
					
				if(j<0)break;
				if(j==0) break;//t0 is zero
				
				Set<Graph> cc = (Set<Graph>) cmap.get(i);
				Set<Graph> tt = (Set<Graph>) tmap.get(j); //kk also contains c
				//System.out.println("We are multiplying c"+i+" and t"+j+" to get c"+m);
				if(cc.isEmpty())System.out.println(i+" Q - cc is empty");
				if(tt.isEmpty())System.out.println(i+" Q - tt is empty");
				int q=0;
				for(Graph gc:cc){
					for(Graph gt:tt){
					
					//	System.out.println("Graph #"+q);q++;
						Graph temp = IEGenerator.fourierMul(gc,gt,1);
						temp = IEGenerator.Relabel(temp).copy();
						//System.out.println("With the fbond");
						temp.putEdge((byte)0,(byte)1);
						//temp = IEGenerator.realMul(temp,fbond).copy();
					//	IEGenerator.GraphDetails(temp);
						c.add(temp);
						h.add(temp);
						
					}
				}
				
			}			
			if(c.isEmpty())System.out.println("c set is empty");
			if(h.isEmpty())System.out.println("h set is empty");
			
			//Adding the Graph Set to cmap
			cmap.add(c);
			hmap.add(h);
			
			//Clearing the sets since it has been added
			c = new HashSet();
			h = new HashSet();
			
		}
		
		
	/*	for(int i=0;i<=n;i++){
			Set<Graph> temp = (Set<Graph>) tmap.get(i);
			System.out.println(" T"+i+" has "+temp.size()+" elements");
		}*/
		
	/*	int p=2;
		Set<Graph> temp = (Set<Graph>) cmap.get(p);
		ClusterViewer.createView("c"+p, temp);*/
	
			
		return cmap;
		
	}
	public Set<Graph> tryout(){
		
		IEGenerator IEGenerator = new IEGenerator();
		Set<Graph> sample = new GraphList<Graph>();
		Relabel Relabel = new Relabel();
		byte[] list = {(byte) 0, (byte) 2, (byte) 1};
		RelabelParameters params = new RelabelParameters(list);
		
		ArrayList chumma = new ArrayList();
		
		
		
		Graph a = GraphFactory.createGraph((byte)3, BitmapFactory.createBitmap((byte)3,false));
		a.coefficient().setNumerator(-2);
		a.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		a.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
		a.putEdge((byte)0, (byte)2);
		a.putEdge((byte)1, (byte)2);
				
		Graph b = GraphFactory.createGraph((byte)2, BitmapFactory.createBitmap((byte)2,true));
		b.coefficient().setNumerator(3);
		b.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		b.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
		
		Graph newGraph = IEGenerator.fourierMul(b,a,1);
			
		System.out.println(newGraph.coefficient());
		Graph y = IEGenerator.Relabel(newGraph).copy();
	
		sample.add(a);
		sample.add(b);
		sample.add(newGraph);
		//sample.add(y);
		
		ClusterViewer.createView("sample",sample);
		
		return null;
		
	}
	
	public static void main(String[] args){
		
		IEGenerator IEGenerator = new IEGenerator();
		
		int n=4;
		
		Graph fbond = GraphFactory.createGraph((byte)2, BitmapFactory.createBitmap((byte)2,true));
		fbond.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		fbond.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
		
		
		Graph a = GraphFactory.createGraph((byte)4, BitmapFactory.createBitmap((byte)4,false));
		a.putEdge((byte)0,(byte)1);
		a.putEdge((byte)1,(byte)2);
		a.putEdge((byte)3,(byte)2);
		a.putEdge((byte)0,(byte)3);
		a.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		a.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
		
		Graph b = GraphFactory.createGraph((byte)3, BitmapFactory.createBitmap((byte)3,false));
		b.putEdge((byte)0,(byte)2);
		b.putEdge((byte)1,(byte)2);
		b.getNode((byte) 0).setType(Metadata.TYPE_NODE_ROOT);
		b.getNode((byte) 1).setType(Metadata.TYPE_NODE_ROOT);
		
		Graph newGraph = IEGenerator.SuperImposeGraphs(a, b);
		Set<Graph> view = new GraphList<Graph>();
		view.add(a);
		view.add(b);
		view.add(newGraph);
		
		
		
		//ClusterViewer.createView("result", view);
		
		//IEGenerator.PYGenerate(n);
		
		IEGenerator.tryout();
		
	}
}
		
/*
 * public Graph ChangeRootNodes(Graph a, int inode1,int inode2,int fnode1,int fnode2){
	
		int newNodeCount=0,node1=0,nodenode2=0;
		newNodeCount = a.nodeCount();
		Graph b = GraphFactory.createGraph((byte) newNodeCount, BitmapFactory.createBitmap((byte)newNodeCount,false));
		IEGenerator IEGenerator = new IEGenerator();
		
		int[] oldlist = new int[a.nodeCount()];
		int[] newlist = new int[a.nodeCount()];
		
		int list = a.nodeCount();
		int iflag=0;
		
		if(inode1==-1){
			for(int i=0;i<a.nodeCount();i++){
				if(a.getNode((byte) i).getType()==Metadata.TYPE_NODE_ROOT && iflag==0 ){
					inode1 = i;
					iflag=1;
				}
				else if(a.getNode((byte) i).getType()==Metadata.TYPE_NODE_ROOT && iflag==1 ){
					inode2 = i;
					break;
				}
		
			}
		}
		System.out.println("old Node 1  = "+inode1+"  ; node 2 = "+inode2);
		
		//Creating old list
		for(int i=0;i<a.nodeCount();i++){
			oldlist[i]=i;
		}
		//Creating new list
		int Duplicationflag=0,dup=0;
		if( inode1==fnode1) { dup =inode1; Duplicationflag =1;}
		else if( inode1 ==fnode2){ dup =inode1; Duplicationflag =1;}
		else if(inode2==fnode1){ dup =inode2; Duplicationflag =1;}
		else if(inode2 ==fnode2 ){ dup =inode2; Duplicationflag =1;}
			
		if(Duplicationflag==1)System.out.println("Duplicate"); 
		else System.out.println("No Duplicate"); 
		
		if(Duplicationflag==0){
			
			for(int i=0;i<oldlist.length;i++){
			
				if(oldlist[i]==inode1)newlist[i]=fnode1;
				else if (oldlist[i]==inode2)newlist[i]=fnode2;
				else if (oldlist[i]==fnode1)newlist[i]=inode1;
				else if (oldlist[i]==fnode2)newlist[i]=inode2;
				else newlist[i]=oldlist[i];
			
			}
		}
		else {
			
			for(int i=0;i<oldlist.length;i++){
				
				if(oldlist[i]==inode1)newlist[i]=fnode1;
				else if (oldlist[i]==inode2)newlist[i]=fnode2;
				else if (oldlist[i]==fnode1){
					if(dup==inode1)newlist[i]=inode2;
					else newlist[i]=inode1;
				}
				else if (oldlist[i]==fnode2){
					if(dup==inode1)newlist[i]=inode2;
					else newlist[i]=inode1;
				}
				else newlist[i]=oldlist[i];
			
			}
			
		}
		
		for(int i=0;i<oldlist.length;i++)
			System.out.println(" old = "+oldlist[i]+" ; Newlist = "+newlist[i]);
		
		//Creating new Graph
		for(int i=0;i<oldlist.length;i++){
			for(int j=0;j<oldlist.length;j++){
				if(a.hasEdge((byte)oldlist[i], (byte)oldlist[j])) {
					b.putEdge((byte)newlist[i], (byte)newlist[j]);
					System.out.println("Adding edges to new Graph between "+newlist[i]+" and "+newlist[j]);
				}
			}
		}
		//Making the required nodes root nodes
		b.getNode((byte) fnode1).setType(Metadata.TYPE_NODE_ROOT);
		b.getNode((byte) fnode2).setType(Metadata.TYPE_NODE_ROOT);
		
		System.out.println(" root end = "+fnode1+" ; rootnode = "+fnode2);
		
		//IEGenerator.GraphDetails(b);
		return b;
	}
 */


