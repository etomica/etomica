package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.operations.DropOrphanNodes;

public class FactorFFT {
	
	int[] factor;
	int[] list;
	IsArticulationPoint p = new IsArticulationPoint();
	
	public void list(int[] array,int time){
		for(int i=0;i<time;i++){
	//		System.out.println(array[i]);
		}
		
	}
	
	public void sort(int[]array,int time){
		
		for (int i=0;i<time;i++)
			for(int j=0;j<time;j++)
				if(array[j+1]<array[j]){
					int t = array[j+1];
					array[j+1] = array[j];
					array[j]=t;
				}
	}
	
	public boolean isArticulated(Graph c, int x, int y){
		
		DropOrphanNodes dropper = new DropOrphanNodes();
		
		Graph g=c.copy();
		
		g.deleteEdge((byte) x,(byte) y);
		
		g=dropper.apply(g);
		
		/*for(int i = 0;(i<g.nodeCount()); i++){       // To find all bonds
			for(int j=0;j<g.getOutDegree( (byte) i);j++){
			System.out.println(i+","+g.getOutNode( (byte) i,(byte) j ));
			}
		}*/
		
		int flag=0;
		//System.out.println("Hello");
		for(int i=0;(i<g.nodeCount()&&(flag==0));i++){
			
			if(p.isArticulationPoint(g,i)){
	//			System.out.println(i+" is an articulation point in factor(isArticulated())");
				flag=1;
			}
	//		if(flag==0)System.out.println(i+" is not an articulation point in factor");
			
		}
		
		if(flag!=0) return true;
		else return false;
		
	}
	
	public boolean isEmpty(Graph g){
		
		Graph c=g.copy();
		int flag=0;
		for(int i=0;((i<c.nodeCount())&&(flag==0));i++)
			for(int j=i+1;((j<c.nodeCount()&&(flag==0)));j++){
				if(c.hasEdge( (byte) i, (byte) j))flag=1;
			}
		
		if (flag==0) return true; // no edges
		else return false;
	}

	/*public boolean factor(Graph g, byte x, byte y){
		
		int factorCounter=0,time=2;
		Graph factor1 = g.copy();
		Graph factor2 = factor1.copy();
		Graph factor3 = factor1.copy();
		int NoArticulationPointFlag=0;
				
		// Deleting bond connecting x and y
		factor1.deleteEdge(x,y);
		
	//	System.out.println(x+","+y+"--->AP");
		
		// Remove all v shaped sub diagrams.
		for(int k=0;k<factor1.nodeCount();k++){  
			if((k==x)||(k==y)) continue;
			else{
			//	System.out.println(k+" = k  ;");
				if((factor1.getOutDegree((byte) k)==2)){
			//		System.out.println(k+" is a point with 2 bonds");
					if ( ((factor1.hasEdge((byte) k,(byte) x))&&(factor1.hasEdge((byte) k, (byte)y))) || ( (factor1.hasEdge((byte) k, (byte)y))&& (factor1.hasEdge((byte) k,(byte) x)))){
			//			System.out.println(k+" has "+factor1.getOutDegree((byte) k)+" bonds");
			//			System.out.println(k+" is a point connected to the AP");
						time++;
							for(int j=factor1.getOutDegree((byte) k)-1;j>=0;j--){
	//							System.out.println("(V) Removing bonds of "+k+"  ("+k+","+factor1.getOutNode((byte) k,(byte)  j)+")");
								factor1.deleteEdge((byte) k,factor1.getOutNode((byte) k,(byte)  j));
								//System.out.println("R");
							}
					}
				}
			}
		}
		
		if(isEmpty(factor1)){	
	//		System.out.println("factor1 is null");
			return true;
		}
				
	/*	if(isDisconnected(factor1));  // for the ring diagrams.
		else {
			System.out.println("This is one big chunk");
			int flag=0;
			for(int j=0;((j<factor1.nodeCount())&&(flag==0));j++){
				if(isArticulationPoint(factor1,j)){
					System.out.println("This chunk has an APT");
					flag=1;
				}
			}
			if(flag!=0) return true;
		}*/
		
	//	while(!isEmpty(factor1)){
			
		/*	for(int i = 0;(i<g.nodeCount()); i++){       // To find all bonds
				for(int j=0;j<g.getOutDegree( (byte) i);j++){
//				System.out.println(i+","+g.getOutNode( (byte) i,(byte) j ));
				}
			}*/
	/*		
			factor2 = factor1.copy();
			factor = new int[10];
			factorCounter=0;
			if(x>y){  // always x needs to be the lower number of the articulation pair
				byte t=x;
				y=t;
				x=y;
			}
			
			int check=0;
			for(int i=0;((i<factor1.nodeCount())&&(check==0));i++){ // get the first node
				if( (factor1.getOutDegree((byte)i)!=0) && (i!=x) && (i!=y) ){  // it should not be any of the AP
					check=1;
	//				System.out.println("The first node of factor is "+i);
					factor[factorCounter]=i;
					factorCounter++;
					time++;
				}
			}
			
		//	System.out.println("factor[0] = "+factor[0]);
		//	System.out.println("factor1.getOutDegree("+factor[0]+") = "+factor1.getOutDegree((byte)factor[0]));
			
			 // add all the connections of the first node
			for(int j= factor1.getOutDegree((byte)factor[0]) - 1;j>=0;j--){
				if ((factor1.getOutNode((byte) factor[0], (byte) j)!=y) && (factor1.getOutNode((byte) factor[0], (byte) j)!=x)){
	//				System.out.println("we add "+factor1.getOutNode((byte) factor[0], (byte) j)+" to the list");
					factor[factorCounter] = factor1.getOutNode((byte) factor[0], (byte) j);
					factorCounter++;
					time++;
				}
			}
			
		//	System.out.println("factor after first's connections is ");
		//	list(factor,factorCounter);
			
			sort(factor,factorCounter-1);
			
		//	System.out.println("factor after first's connections sorting is ");
		//	list(factor,factorCounter);
			
			
			// Getting one all the connections of the nodes connected to the first node and so on..
			for(int i=0;i<factorCounter+1;i++){  
				
				if((i==x)||(i==y)){
		//			System.out.println("AP encountered so leave this iteration..");	
					continue;
				}
				
		//		System.out.println(i+" = i and factorCounter ="+factorCounter);
				
				int a = factor[i],flag=0;
			
		//		System.out.println("Checking connections with "+a);
			
				for(int j=factor1.getOutDegree((byte) a)-1;j>=0;j--){
				
				//	System.out.println(a+" has "+ factor1.getOutDegree((byte) a)+" bonds");
				
					for(int p=0;p<factorCounter;p++){
						if(factor[p]==factor1.getOutNode((byte)a,(byte)j)){
							flag=1; 
		//					System.out.println(factor1.getOutNode((byte)a,(byte)j)+" is already in the factor so DONT ADD");
						}
					}
				
					if((flag==0)&&(x!=factor1.getOutNode((byte)a,(byte)j))&&(y!=factor1.getOutNode((byte)a,(byte)j))){
		//				System.out.println(factor1.getOutNode((byte)a,(byte)j)+" is not in the factor and not one of the AP so PLEASE ADD");
						factor[factorCounter]=factor1.getOutNode((byte)a,(byte)j);
						factorCounter++;
						time++;
						sort(factor,factorCounter-1);
						i=0;		// since we need to start from beginning
					}
					
					flag=0;
				
				}	
			}
			
			// add the articulation pair to the list
			factor[factorCounter]=x;factorCounter++;
			factor[factorCounter]=y;factorCounter++;
			
			sort(factor,factorCounter-1);
			
	//		System.out.println("factor after finding the complete connected list is ");
			for(int i=0;i<factorCounter;i++){
	//			System.out.println(factor[i]);
			}
			
			for(int i=0;i<factor1.nodeCount();i++){  // deleting node edges for nodes not in factor1
			
				int delflag=0;
				for(int j=0;((j<factorCounter)&&(delflag==0));j++){
					if(i==factor[j]){
	//					System.out.println(i+" is in the list so dont delete");
						delflag=1;
					}
				}
				
				if((delflag==0)&&(i!=x)&&(i!=y)){
	//				System.out.println("Deleting "+i+" since its not in factor");
					for(int k=factor1.getOutDegree((byte) i)-1;k>=0;k--){
						factor1.deleteEdge((byte) i,factor1.getOutNode((byte) i,(byte)  k));
					}
				}
			}	
			
			if(!isArticulated(factor1)){  // no articulation point in factor
	//			System.out.println("No articulation points in factor1");
				NoArticulationPointFlag=1;
			}
//			else System.out.println("Factor has articulation Point");
			
			factor1=factor2.copy(); 
			
			/*for(int i = 0;(i<g.nodeCount()); i++){       // To find all bonds
				for(int j=0;j<g.getOutDegree( (byte) i);j++){
//				System.out.println(i+","+g.getOutNode( (byte) i,(byte) j ));
				}
			}*/
			
		/*	for(int i=0;i<factorCounter;i++){
				int a = factor[i];
				if((a!=x)&&(a!=y)){
//					System.out.println("Deleting edges of nodes "+a+ " since its not an AP");
					for(int k=factor1.getOutDegree((byte) a)-1;k>=0;k--){
						factor1.deleteEdge((byte) a,factor1.getOutNode((byte) a,(byte)  k));
					}
				}
			}
			
		}
		
		if(NoArticulationPointFlag==0)return true;
		else return false;
	}	*/
	
public boolean factor(Graph g, byte x, byte y){
		
		int factorCounter=0;
		Graph factor1 = g.copy();
		Graph factor2 = factor1.copy();
		int NoArticulationPointFlag=0;
				
		while(!isEmpty(factor1)){
			
			factor1.putEdge((byte) x,(byte) y);
			
			factor2 = factor1.copy();
			factor = new int[10];
			factorCounter=0;
		//	System.out.println("x = "+x+";   y = "+y);
			if(x>y){  // always x needs to be the lower number of the articulation pair
				byte t=x;
				x=y;
				y=t;
		//		System.out.println("x = "+x+";   y = "+y);
			}
			
			
			int check=0;
			for(int i=0;((i<factor1.nodeCount())&&(check==0));i++){ // get the first node
				if( (factor1.getOutDegree((byte)i)!=0) && (i!=x) && (i!=y) ){  // it should not be any of the AP
					check=1;
		//			System.out.println("The first node of factor is "+i);
					factor[factorCounter]=i;
					factorCounter++;
					
				}
			}
			
		//	System.out.println("factor[0] = "+factor[0]);
		//	System.out.println("factor1.getOutDegree("+factor[0]+") = "+factor1.getOutDegree((byte)factor[0]));
			
			 // add all the connections of the first node
			for(int j= factor1.getOutDegree((byte)factor[0])-1;j>=0;j--){
				if ((factor1.getOutNode((byte) factor[0], (byte) j)!=y) && (factor1.getOutNode((byte) factor[0], (byte) j)!=x)){
		//			System.out.println("since the number is not AP we add "+factor1.getOutNode((byte) factor[0], (byte) j)+" to the list");
					factor[factorCounter] = factor1.getOutNode((byte) factor[0], (byte) j);
					factorCounter++;
				}
			}
			
		//	System.out.println("factor after first's connections is ");
		//	list(factor,factorCounter);
			
			sort(factor,factorCounter-1);
			
	//		System.out.println("factor after sorting first's connections  is ");
	//		list(factor,factorCounter);
			
			
			// Getting all the connections of the nodes connected to the first node and so on..
			for(int i=0;i<=factorCounter;i++){  
				
				int a = factor[i],flag=0;
				
				if((a==x)||(a==y)){
	//				System.out.println("AP encountered so leave this iteration..");	
					continue;
				}
				
			//	System.out.println(i+" = i and factorCounter ="+factorCounter);
				
			//	System.out.println("Checking connections with "+a);
			
				for(int j=factor1.getOutDegree((byte) a)-1;j>=0;j--){
				
					//System.out.println(a+" has "+ factor1.getOutDegree((byte) a)+" bonds");
				
					for(int p=0;p<factorCounter;p++){
						if(factor[p]==factor1.getOutNode((byte)a,(byte)j)){
							flag=1; 
	//						System.out.println(factor1.getOutNode((byte)a,(byte)j)+" is already in the factor so DONT ADD");
						}
					}
				
					if((flag==0)&&(x!=factor1.getOutNode((byte)a,(byte)j))&&(y!=factor1.getOutNode((byte)a,(byte)j))){
	//					System.out.println(factor1.getOutNode((byte)a,(byte)j)+" is not in the factor and not one of the AP so PLEASE ADD");
						factor[factorCounter]=factor1.getOutNode((byte)a,(byte)j);
						factorCounter++;
						sort(factor,factorCounter-1);
						i=0;		// since we need to start from beginning
					}
					
					flag=0;
				
				}	
			}
			
			// add the articulation pair to the list
			factor[factorCounter]=x;factorCounter++;
			factor[factorCounter]=y;factorCounter++;
			
			sort(factor,factorCounter-1);
			
	//		System.out.println("factor after finding the complete connected list is ");
	/*		for(int i=0;i<factorCounter;i++){
				System.out.println(factor[i]);
			}*/
			
			for(int i=0;i<factor1.nodeCount();i++){  // deleting node edges for nodes not in factor1
			
				int delflag=0;
				for(int j=0;((j<factorCounter)&&(delflag==0));j++){
					if(i==factor[j]){
	//					System.out.println(i+" is in the list so dont delete");
						delflag=1;
					}
				}
				
				if((delflag==0)&&(i!=x)&&(i!=y)){
	//				System.out.println("Deleting "+i+" since its not in factor");
					for(int k=factor1.getOutDegree((byte) i)-1;k>=0;k--){
						factor1.deleteEdge((byte) i,factor1.getOutNode((byte) i,(byte)  k));
					}
				}
			}	
			
			if(!isArticulated(factor1,x,y)){  // no articulation point in factor
	//			System.out.println("No articulation points in factor1");
				NoArticulationPointFlag=1;
				return false;
			}
	//		else System.out.println("Factor has articulation Point");
			
			factor1=factor2.copy(); 
			
			for(int i=0;i<factorCounter;i++){
				int a = factor[i];
				if((a!=x)&&(a!=y)){
	//				System.out.println("Deleting nodes present in factor like "+a+ " since its not an AP");
					for(int k=factor1.getOutDegree((byte) a)-1;k>=0;k--){
						factor1.deleteEdge((byte) a,factor1.getOutNode((byte) a,(byte)  k));
					}
				}
			}
			
			factor1.deleteEdge((byte) x,(byte) y);
			
		}
		
		if(NoArticulationPointFlag==0)return true;
		else return false;
	}	
}
