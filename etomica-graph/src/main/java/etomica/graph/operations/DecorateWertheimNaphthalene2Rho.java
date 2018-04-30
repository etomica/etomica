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
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;

/**
 * modified from Hye Min's DecorateWertheim2SiteRho class
 * Performs decoration of diagrams by various density points.
 * 2 site model has  rho0, rhoA, rhoB and rhoAB.
 * naphthalene Wertheim Theory
 * @author shu
 * Date:Sep.9.2011
 * 
 */
public class DecorateWertheimNaphthalene2Rho implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
      char capAABond = ((DecorateNaphthalene2Site)params).capAABond;
      char capABBond = ((DecorateNaphthalene2Site)params).capABBond;
      char capBABond = ((DecorateNaphthalene2Site)params).capBABond;
      char capBBBond = ((DecorateNaphthalene2Site)params).capBBBond;

      Set<Graph> rhoA = ((DecorateNaphthalene2Site)params).rhoA;
      Set<Graph> rhoB = ((DecorateNaphthalene2Site)params).rhoB;
      Set<Graph> rhoAB = ((DecorateNaphthalene2Site)params).rhoAB;
      MulFlexibleParameters mfp = ((DecorateNaphthalene2Site)params).mfp;
      MulScalar mulScalar = new MulScalar();
      Set<Graph> rho0CA = new HashSet<Graph>();
      Set<Graph> rho0CB = new HashSet<Graph>();
      Set<Graph> rho0CAB = new HashSet<Graph>();
      MulScalarParameters msp = new MulScalarParameters(-1,1);
      MulFlexible mulFlex = new MulFlexible();
      rho0CA.addAll(mulScalar.apply(rhoA, msp));
      rho0CB.addAll(mulScalar.apply(rhoB, msp));
      rho0CAB.addAll(mulScalar.apply(rhoA, msp));
      rho0CAB.addAll(mulScalar.apply(rhoB, msp));
      rho0CAB.addAll(mulScalar.apply(rhoAB, msp));
      
      char sigmaAB = 'a';
      char sigmaA = 'b';
      char sigmaB = 'c';
      char sigma0 = 'A';
      char sigmaABFinal = 'X';
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
      Graph sigmaAPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single point
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
      Graph sigmaBPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single point
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
      
      while (!pool.isEmpty()){

    	  Graph g = pool.remove(pool.size()-1).copy();//remove last element
        
        if(g.factors().length == 0){
              g.setNumFactors(4);
              g.addFactors(new int[]{g.nodeCount(),0,0,0});
        }
          
        
        // check the type of bond that links the point; loop over all points
        
    	int foundSigmaPoint = 0;
    	int[] bitMap = new int[4];
    	bitMap[0] = maskAB;
    	bitMap[maskA] = maskB;
    	bitMap[maskB] = maskA;
        for (Node node1 : g.nodes()) {
        	int bits1 = 0;
        	for (Node node2 : g.nodes()) {
        		if (node1.getId() == node2.getId())continue;
        		if (!g.hasEdge(node1.getId(), node2.getId()))continue;
	              Edge edge12 = g.getEdge(node1.getId(),node2.getId());
	              if (edge12.getColor() == capAABond || edge12.getColor() == capABBond){
	            	  bits1 = (bits1 | maskA);//bits1=bits+maskA(4)
	              } else if (edge12.getColor() == capBABond || edge12.getColor() == capBBBond){
	            	  bits1  = (bits1  | maskB);//bits1=bits+maskB(2)
		            	
	              }		
        	}
        	char color1 = node1.getColor();
        	if (node1.getType() != Metadata.TYPE_NODE_ROOT && (foundSigmaPoint & bitMap[bits1])== 0 && (color1 ==sigma0)){
        		g.getNode(node1.getId()).setColor(sigma[bits1]);
        		foundSigmaPoint|= bitMap[bits1];
     
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

  
  //inner class
  public static class DecorateNaphthalene2Site implements Parameters {
	  public final MulFlexibleParameters mfp;
	  public final char capAABond, capABBond,capBABond, capBBBond;
	  public final Set<Graph> rhoA,rhoB, rhoAB;
	  public DecorateNaphthalene2Site(MulFlexibleParameters mfp, char capAABond, char capABBond,char capBABond, char capBBBond, Set<Graph> rhoA, Set<Graph> rhoB, Set<Graph> rhoAB) {
		  this.capAABond = capAABond;
		  this.capABBond = capABBond;
		  this.capBABond = capBABond;
		  this.capBBBond = capBBBond;
		  
		  this.rhoA = rhoA;
		  this.rhoB = rhoB;
		  this.rhoAB = rhoAB;
		  this.mfp = mfp;
	  }
	  
  }
}