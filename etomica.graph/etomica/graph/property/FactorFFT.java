/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.property;

import etomica.graph.model.Graph;
import etomica.graph.operations.DropOrphanNodes;

public class FactorFFT {
	
	IsArticulationPoint p = new IsArticulationPoint();
	
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

		for(int i=0;i<g.nodeCount();i++){
			
			if(p.isArticulationPoint(g,i)){
				return true;
			}
			
		}
		
		return false;
		
	}
	
	public boolean isEmpty(Graph g){
		return g.edgeCount() == 0;
	}
	public boolean factor(Graph g, byte x, byte y){
		
		Graph factor1 = g.copy();
    if(x>y){  // always x needs to be the lower number of the articulation pair
      byte t=x;
      x=y;
      y=t;
    }
				
		while(!isEmpty(factor1)){
			
			factor1.putEdge(x,y);
			
			Graph factor2 = factor1.copy();
			int[] factor = new int[10];
			int factorCounter=0;
			
			for(int i=0;i<factor1.nodeCount();i++){ // get the first node
				if( (factor1.getOutDegree((byte)i)!=0) && (i!=x) && (i!=y) ){  // it should not be any of the AP
		//			System.out.println("The first node of factor is "+i);
					factor[factorCounter]=i;
					factorCounter++;
					break;
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
if (true) throw new RuntimeException("articulation pair shouldn't be in here");
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
			
			factor1.deleteEdge(x, y);
			
		}
		
		return true;
	}	
}
