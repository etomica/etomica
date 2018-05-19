/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graph.operations;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;

import etomica.graph.model.BitmapFactory;
import etomica.graph.model.Edge;
import etomica.graph.model.Graph;
import etomica.graph.model.GraphFactory;
import etomica.graph.model.GraphList;
import etomica.graph.model.Metadata;
import etomica.graph.model.Node;
import etomica.graph.model.comparators.ComparatorBiConnected;
import etomica.graph.model.comparators.ComparatorChain;
import etomica.graph.model.comparators.ComparatorNumEdges;
import etomica.graph.model.comparators.ComparatorNumFieldNodes;
import etomica.graph.model.comparators.ComparatorNumNodes;
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;

/**
 * Performs decoration of diagrams by various density points.
 * 2 site model has  rho0, rhoA, rhoB and rhoAB.
 *
 * @author Hye Min Kim
 */
public class DecorateWertheim2Site implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
      int n = ((DecorateWertheimParameters2Site)params).n;
      char capFABBond = ((DecorateWertheimParameters2Site)params).capFABBond;
      char capFBABond = ((DecorateWertheimParameters2Site)params).capFBABond;
      Set<Graph> cA = ((DecorateWertheimParameters2Site)params).cA;
      Set<Graph> cB = ((DecorateWertheimParameters2Site)params).cB;
      Set<Graph> cAB = ((DecorateWertheimParameters2Site)params).cAB;
      MulFlexibleParameters mfp = ((DecorateWertheimParameters2Site)params).mfp;
      MulScalar mulScalar = new MulScalar();
      Set<Graph> rho0CA = new HashSet<Graph>();
      Set<Graph> rho0CB = new HashSet<Graph>();
      Set<Graph> rho0CAB = new HashSet<Graph>();
      MulScalarParameters msp = new MulScalarParameters(-1,1);
      MulFlexible mulFlex = new MulFlexible();
      char[] flexColors = new char[0];
      MulFlexibleParameters mfpCi = MulFlexibleParameters.makeParametersOnlyRootPt(flexColors, (byte)(n-1));
      IsoFree isoFree = new IsoFree();
      Set<Graph> cAcB = mulFlex.apply(cA, cB, mfpCi);
      cAcB = isoFree.apply(cAcB,null);
      rho0CA.addAll(mulScalar.apply(cA, msp));
      rho0CB.addAll(mulScalar.apply(cB, msp));
      rho0CAB.addAll(mulScalar.apply(cA, msp));
      rho0CAB.addAll(mulScalar.apply(cB, msp));
      rho0CAB.addAll(mulScalar.apply(cAcB, msp));
      rho0CAB.addAll(mulScalar.apply(cAB, msp));
      ComparatorChain comp = new ComparatorChain();
      comp.addComparator(new ComparatorNumFieldNodes());
      comp.addComparator(new ComparatorBiConnected());
      comp.addComparator(new ComparatorNumEdges());
      comp.addComparator(new ComparatorNumNodes());
      GraphList graphList = new GraphList(comp);
      graphList.addAll(rho0CAB);
      
  	char sigmaAB = 'a';//color of node
  	char sigmaA = 'b';
  	char sigmaB = 'c';
  	char sigma0 = 'A';
  	char sigmaABFinal = 'X';//color of node
  	char sigmaAFinal = 'B';
  	char sigmaBFinal = 'C';
  	char sigma0Final = 'D';
  	char[] sigma = new char[]{sigmaAB,sigmaA,sigmaB,sigma0};
  	int maskAB = 0x4;//4
  	int maskA = 0x2;//2
  	int maskB = 0x1;//1

  
      for (Graph g:rho0CA){
          g.setNumFactors(4);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaA);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0});
      }
      Graph sigmaAPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaAPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaAPoint.getNode((byte)0).setColor(sigmaA);
      sigmaAPoint.setNumFactors(4);
      sigmaAPoint.addFactors(new int[]{0,1,0,0});//1 sigmaA
      rho0CA.add(sigmaAPoint);
      
      for (Graph g:rho0CB){
          g.setNumFactors(4);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaB);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0});
      }
      Graph sigmaBPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaBPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaBPoint.getNode((byte)0).setColor(sigmaB);
      sigmaBPoint.setNumFactors(4);
      sigmaBPoint.addFactors(new int[]{0,0,1,0});//1 sigmaB
      rho0CB.add(sigmaBPoint);
     
      for (Graph g:rho0CAB){
          g.setNumFactors(4);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaAB);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0});
      }
     
      Graph sigmaABPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaABPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaABPoint.getNode((byte)0).setColor(sigmaAB);
      sigmaABPoint.setNumFactors(4);
      sigmaABPoint.addFactors(new int[]{0,0,0,1});//1 sigmaAB
      rho0CAB.add(sigmaABPoint);
       
      Set<Graph> result = new HashSet<Graph>();    
      ArrayList<Graph> pool = new ArrayList<Graph>();
      Set<Graph> product = new GraphList();
      pool.addAll(argument);
      CoefficientImpl coefficient = new CoefficientImpl(0);
      
      while (!pool.isEmpty()){//argument = p
        Graph g0 = pool.get(pool.size()-1);
        if (argument.contains(g0)&&(g0.edgeCount()==g0.nodeCount()-1)){
        	System.out.println("old coefficient "+coefficient);
        	coefficient = new CoefficientImpl(0);
        	System.out.println("g0 "+g0);
        }
    	  Graph g = pool.remove(pool.size()-1).copy();//remove last element
        
        if(g.factors().length == 0){
              g.setNumFactors(4);//sigma0, sigmaA, sigmaB, sigmaAB
              g.addFactors(new int[]{g.nodeCount(),0,0,0});
        }
          
        int[] bits = new int[g.nodeCount()];
    	int foundSigmaPoint = 0;
    	int counterNodeDecorate = 0;
    	int[] bitMap = new int[4];
    	bitMap[0] = maskAB;
    	bitMap[maskA] = maskB;
    	bitMap[maskB] = maskA;
        for (Node node1 : g.nodes()) {
        	int id1 = node1.getId();
        	for (Node node2 : g.nodes()) {
        		int id2 = node2.getId();
        		if (node1.getId() >= node2.getId())continue;
        		if (!g.hasEdge(node1.getId(), node2.getId()))continue;
	              Edge edge12 = g.getEdge(node1.getId(),node2.getId());
	              if (edge12.getColor() == capFABBond){
	            	  bits[id1] = (bits[id1]  | maskA);//bits1=bits+maskA(4)
	              } else if (edge12.getColor() == capFBABond){
	            	  bits[id1]  = (bits[id1]  | maskB);//bits1=bits+maskB(2)
	
	              }
	              if (edge12.getColor() == capFBABond){
	            	  bits[id2]  = (bits[id2] | maskA);
	            	  
	              } else if (edge12.getColor() == capFABBond){
	            	  bits[id2]  = (bits[id2] | maskB);
	            	  
	              }		
        	}
        	char color1 = node1.getColor();
        	if (node1.getType() != Metadata.TYPE_NODE_ROOT && (foundSigmaPoint & bitMap[bits[id1]])== 0 && (color1 ==sigma0)){
        		g.getNode(node1.getId()).setColor(sigma[bits[id1]]);
        		foundSigmaPoint|= bitMap[bits[id1]];
        		if(sigma[bits[id1]] != sigma0){
        			counterNodeDecorate++;
        		}
        	}
        }
        int sigma0Factor = 0;
        int sigmaAFactor = 0;
        int sigmaBFactor = 0;
        int sigmaABFactor = 0;

        for(Node node1 : g.nodes()){
      	  if (node1.getType() == Metadata.TYPE_NODE_ROOT){
    		  node1.setColor('E');
    	  }
      	   else if(node1.getColor()==sigmaABFinal){
        		sigmaABFactor++;
        	}
        	else if(node1.getColor()==sigmaAFinal){
        		sigmaAFactor++;
        	}
        	else if(node1.getColor()==sigmaBFinal){
        		sigmaBFactor++;
        	}
        	else if(node1.getColor()==sigma0Final||node1.getColor()==sigma0){
        		sigma0Factor++;
        	}
        }

          if(foundSigmaPoint == 0){
        	  if(g.factors().length != 4) {
        		  g.setNumFactors(4);
        		  g.addFactors(new int[]{g.nodeCount(),0,0,0});
        	  }
              result.add(g);
              continue;
          }
          
          g.setNumFactors(4);
          g.addFactors(new int[]{sigma0Factor,sigmaAFactor,sigmaBFactor,sigmaABFactor});//decorate only 1 point at one time
          product.clear();
          product.add(g);
          if ((foundSigmaPoint & maskA) != 0){
	          product = mulFlex.apply(product, rho0CA, mfp);
	          for (Graph gP:product){
	              for(Node node1 : gP.nodes()){
	                  if(node1.getColor() == sigmaA){
	                      if(gP.nodeCount()== g.nodeCount())
	                          node1.setColor(sigmaAFinal);//sigmaA
	                      else node1.setColor(sigma0);//sigma0
	                  }
	              }
	          }
          }
          if ((foundSigmaPoint & maskB) != 0){
        	  Set<Graph> g1 = new HashSet<Graph>();//single diagram
              Set<Graph> product1 = new GraphList();
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CB, mfp);
    	          for (Graph gP2:newProduct){
    	              for(Node node1 : gP2.nodes()){
    	                  if(node1.getColor() == sigmaB){
    	                      if(gP2.nodeCount()== gP.nodeCount())
    	                          node1.setColor(sigmaBFinal);//sigmaB
    	                      else node1.setColor(sigma0);//sigma0
    	                  }
    	              }
    	          }
    	          product1.addAll(newProduct);
              }
              product = product1;

          }
          if ((foundSigmaPoint & maskAB) != 0){
        	  Set<Graph> g1 = new HashSet<Graph>();//single diagram
              Set<Graph> product1 = new GraphList();
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CAB, mfp);
    	          for (Graph gP2:newProduct){
    	              for(Node node1 : gP2.nodes()){
    	                  if(node1.getColor() == sigmaAB){
    	                      if(gP2.nodeCount()== gP.nodeCount())
    	                          node1.setColor(sigmaABFinal);//sigmaAB
    	                      else node1.setColor(sigma0);//sigma0
    	                  }
    	              }
    	          }
    	          product1.addAll(newProduct);
              }
              product = product1;
          }

          pool.addAll(product);
          for (Graph gP:product){
        	  if (gP.nodeCount() == gP.edgeCount()+1 && !graphHasEdgeColor(gP, 'B')){
        		  if (gP.nodeCount() == 4){
        			  boolean branch = false;
        			  boolean sameColor = false;
        		  for (Node node: gP.nodes()){
        			  if (gP.getOutDegree(node.getId())==3){
        				  branch = true;
        			  } else if(gP.getOutDegree(node.getId()) == 2){
        				  byte id0 = gP.getOutNode(node.getId(), (byte)0);
        				  byte id1 = gP.getOutNode(node.getId(), (byte)1);
        				  //boolean sameDirection = (id0-node.getId())*(id1-node.getId())>0;//where the given bond is going to.. index = number of bond
        				  char color0 = gP.getEdge(node.getId(), id0).getColor();
        				  char color1 = gP.getEdge(node.getId(), id1).getColor();
        				  sameColor = color0 == color1;
        				  if(sameColor){
        					  break;
        				  }
        				  
        				  
        			  }
        		  }
        		  if(!branch&&!sameColor){
        			  if (g.nodeCount() == gP.nodeCount()){
        				  coefficient.add(gP.coefficient());
        			  }
        			  //System.out.println("g "+g);
        			  //System.out.println("gP "+gP);
        		  }

        		  }
                  //break;
        	  }
          }

      }
      for (Graph g:result){
          for(Node node1 : g.nodes()){
              node1.setColor('A');
          }
      }
      return result;
  }
  
  public static boolean graphHasEdgeColor(Graph g, char color) {
      for (Edge edge : g.edges()) {
          if (edge.getColor() == color) {
              return true;
          }
      }
      return false;
  }
  
  public static int countEdgeColor(Graph g, char color){
	  int count = 0;
	  for(Edge edge:g.edges()){
		  if(edge.getColor() == color){
			  count++;
		  }
	  }
	  return count;
  }

  public static class DecorateWertheimParameters2Site implements Parameters {
    public final MulFlexibleParameters mfp;
    public final char capFABBond, capFBABond;
    public final Set<Graph> cA,cB, cAB;
    public final int n;
    public DecorateWertheimParameters2Site(int n, MulFlexibleParameters mfp, char capFABBond, char capFBABond, Set<Graph> cA, Set<Graph> cB, Set<Graph> cAB) {
      this.n=n;   
      this.capFABBond = capFABBond;
      this.capFBABond = capFBABond;
      this.cA = cA;
      this.cB = cB;
      this.cAB = cAB;
      this.mfp = mfp;
    }
   
  }
}