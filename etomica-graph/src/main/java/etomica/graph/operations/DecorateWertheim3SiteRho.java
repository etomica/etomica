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
import etomica.graph.model.impl.CoefficientImpl;
import etomica.graph.operations.MulFlexible.MulFlexibleParameters;
import etomica.graph.property.NumFieldNodes;

/**
 * Performs decoration of diagrams by various density points.
 * 3 site model has  rho0, rhoA, rhoB, rhoC, rhoAB, rhoAC, rhoBC and rhoABC.
 *
 * @author Hye Min Kim
 */
public class DecorateWertheim3SiteRho implements Unary {

  public Set<Graph> apply(Set<Graph> argument, Parameters params) {
      char capFACBond = ((DecorateWertheimParameters3Site)params).capFACBond;
      char capFBCBond = ((DecorateWertheimParameters3Site)params).capFBCBond;
      char capFCABond = ((DecorateWertheimParameters3Site)params).capFCABond;
      char capFCBBond = ((DecorateWertheimParameters3Site)params).capFCBBond;
      char mCapFACBond = ((DecorateWertheimParameters3Site)params).mCapFACBond;
      char mCapFBCBond = ((DecorateWertheimParameters3Site)params).mCapFBCBond;
      char mCapFCABond = ((DecorateWertheimParameters3Site)params).mCapFCABond;
      char mCapFCBBond = ((DecorateWertheimParameters3Site)params).mCapFCBBond;
      Set<Graph> rhoA = ((DecorateWertheimParameters3Site)params).rhoA;
      Set<Graph> rhoB = ((DecorateWertheimParameters3Site)params).rhoB;
      Set<Graph> rhoC = ((DecorateWertheimParameters3Site)params).rhoC;
      Set<Graph> rhoAB = ((DecorateWertheimParameters3Site)params).rhoAB;
      Set<Graph> rhoAC = ((DecorateWertheimParameters3Site)params).rhoAC;
      Set<Graph> rhoBC = ((DecorateWertheimParameters3Site)params).rhoBC;
      Set<Graph> rhoABC = ((DecorateWertheimParameters3Site)params).rhoABC;
      MulFlexibleParameters mfp = ((DecorateWertheimParameters3Site)params).mfp;
      MulScalar mulScalar = new MulScalar();
      Set<Graph> rho0CA = new HashSet<Graph>();
      Set<Graph> rho0CB = new HashSet<Graph>();
      Set<Graph> rho0CC = new HashSet<Graph>();
      Set<Graph> rho0CAB = new HashSet<Graph>();
      Set<Graph> rho0CAC = new HashSet<Graph>();
      Set<Graph> rho0CBC = new HashSet<Graph>();
      Set<Graph> rho0CABC = new HashSet<Graph>();
      MulScalarParameters msp = new MulScalarParameters(-1,1);
      MulFlexible mulFlex = new MulFlexible();
      rho0CA.addAll(mulScalar.apply(rhoA, msp));
      rho0CB.addAll(mulScalar.apply(rhoB, msp));
      rho0CC.addAll(mulScalar.apply(rhoC, msp));
      rho0CAB.addAll(mulScalar.apply(rhoA, msp));
      rho0CAB.addAll(mulScalar.apply(rhoB, msp));
      rho0CAB.addAll(mulScalar.apply(rhoAB, msp));
      rho0CAC.addAll(mulScalar.apply(rhoA, msp));
      rho0CAC.addAll(mulScalar.apply(rhoC, msp));
      rho0CAC.addAll(mulScalar.apply(rhoAC, msp));
      rho0CBC.addAll(mulScalar.apply(rhoB, msp));
      rho0CBC.addAll(mulScalar.apply(rhoC, msp));
      rho0CBC.addAll(mulScalar.apply(rhoBC, msp));     
      rho0CABC.addAll(mulScalar.apply(rhoA, msp));
      rho0CABC.addAll(mulScalar.apply(rhoB, msp));
      rho0CABC.addAll(mulScalar.apply(rhoC, msp));
      rho0CABC.addAll(mulScalar.apply(rhoAB, msp));
      rho0CABC.addAll(mulScalar.apply(rhoAC, msp));
      rho0CABC.addAll(mulScalar.apply(rhoBC, msp));
      rho0CABC.addAll(mulScalar.apply(rhoABC, msp));
      
  	char sigmaABC = 'a';//color of node
  	char sigmaAB = 'b'; 
  	char sigmaAC = 'c';
  	char sigmaBC = 'd';
  	char sigmaA = 'e';
  	char sigmaB = 'f';
  	char sigmaC = 'g';
  	char sigma0 = 'A';
  	char sigmaABCFinal = 'X';
  	char sigmaABFinal = 'B'; 
  	char sigmaACFinal = 'C';
  	char sigmaBCFinal = 'D';
  	char sigmaAFinal = 'E';
  	char sigmaBFinal = 'F';
  	char sigmaCFinal = 'G';
  	char sigma0Final = 'H';
  	char[] sigma = new char[]{sigmaABC,sigmaAB,sigmaAC,sigmaA,sigmaBC,sigmaB,sigmaC,sigma0};
  	int maskA = 0x4;
  	int maskB = 0x2;
  	int maskC = 0x1;
  	int maskAB = 0x8;
  	int maskAC = 0x10;
  	int maskBC = 0x20;
  	int maskABC = 0x40;//64

  
      for (Graph g:rho0CA){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaA);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
      Graph sigmaAPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaAPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaAPoint.getNode((byte)0).setColor(sigmaA);
      sigmaAPoint.setNumFactors(8);
      sigmaAPoint.addFactors(new int[]{0,1,0,0,0,0,0,0});//1 sigmaA
      rho0CA.add(sigmaAPoint);
      
      for (Graph g:rho0CB){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaB);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
      Graph sigmaBPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaBPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaBPoint.getNode((byte)0).setColor(sigmaB);
      sigmaBPoint.setNumFactors(8);
      sigmaBPoint.addFactors(new int[]{0,0,1,0,0,0,0,0});//1 sigmaB
      rho0CB.add(sigmaBPoint);
      
      for (Graph g:rho0CC){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaC);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
     
      Graph sigmaCPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaCPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaCPoint.getNode((byte)0).setColor(sigmaC);
      sigmaCPoint.setNumFactors(8);
      sigmaCPoint.addFactors(new int[]{0,0,0,1,0,0,0,0});//1 sigmaC
      rho0CC.add(sigmaCPoint);
     
      for (Graph g:rho0CAB){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaAB);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
     
      Graph sigmaABPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaABPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaABPoint.getNode((byte)0).setColor(sigmaAB);
      sigmaABPoint.setNumFactors(8);
      sigmaABPoint.addFactors(new int[]{0,0,0,0,1,0,0,0});//1 sigmaAB
      rho0CAB.add(sigmaABPoint);
     
      for (Graph g:rho0CAC){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaAC);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
     
      Graph sigmaACPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaACPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaACPoint.getNode((byte)0).setColor(sigmaAC);
      sigmaACPoint.setNumFactors(8);
      sigmaACPoint.addFactors(new int[]{0,0,0,0,0,1,0,0});//1 sigmaAC
      rho0CAC.add(sigmaACPoint);
      
      for (Graph g:rho0CBC){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaBC);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
     
      Graph sigmaBCPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaBCPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaBCPoint.getNode((byte)0).setColor(sigmaBC);
      sigmaBCPoint.setNumFactors(8);
      sigmaBCPoint.addFactors(new int[]{0,0,0,0,0,0,1,0});//1 sigmaBC
      rho0CBC.add(sigmaBCPoint);
     
      for (Graph g:rho0CABC){
          g.setNumFactors(8);
          for(Node node1 : g.nodes()){
              if(node1.getType()==Metadata.TYPE_NODE_ROOT){
                  node1.setColor(sigmaABC);
              }
          }
         
          g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
      }
      
      Graph sigmaABCPoint = GraphFactory.createGraph((byte)1, BitmapFactory.createBitmap((byte)1,true));//add single pt
      sigmaABCPoint.getNode((byte)0).setType(Metadata.TYPE_NODE_ROOT);
      sigmaABCPoint.getNode((byte)0).setColor(sigmaABC);
      sigmaABCPoint.setNumFactors(8);
      sigmaABCPoint.addFactors(new int[]{0,0,0,0,0,0,0,1});//1 sigmaABC
      rho0CABC.add(sigmaABCPoint);

      Set<Graph>[] rho0CASets = splitSets(rho0CA, mfp.nFieldPoints);
      Set<Graph>[] rho0CBSets = splitSets(rho0CB, mfp.nFieldPoints);
      Set<Graph>[] rho0CCSets = splitSets(rho0CC, mfp.nFieldPoints);
      Set<Graph>[] rho0CABSets = splitSets(rho0CAB, mfp.nFieldPoints);
      Set<Graph>[] rho0CBCSets = splitSets(rho0CBC, mfp.nFieldPoints);
      Set<Graph>[] rho0CACSets = splitSets(rho0CAC, mfp.nFieldPoints);
      Set<Graph>[] rho0CABCSets = splitSets(rho0CABC, mfp.nFieldPoints);

      Set<Graph> result = new HashSet<Graph>();    
      ArrayList<Graph> pool = new ArrayList<Graph>();
      Set<Graph> product = new GraphList(null);
      pool.addAll(argument);
      CoefficientImpl coefficient = new CoefficientImpl(0);
      
      while (!pool.isEmpty()){//argument = p
          Graph g0 = pool.get(pool.size()-1);
          if (argument.contains(g0)&&(g0.edgeCount()==g0.nodeCount()-1)){
          	//System.out.println("old coefficient "+coefficient);
          	coefficient = new CoefficientImpl(0);
          	//System.out.println("g0 "+g0);
          }
        Graph g = pool.remove(pool.size()-1).copy();//remove last element
        
        if(g.factors().length == 0){
              g.setNumFactors(8);//sigma0, sigmaA, sigmaB, sigmaC, sigmaAB, sigmaAC, sigmaBC, sigmaABC
              g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
        }
          
        int[] bits = new int[g.nodeCount()];
    	int foundSigmaPoint = 0;
    	int[] bitMap = new int[8];
    	bitMap[0] = maskABC;
    	bitMap[maskA] = maskBC;
    	bitMap[maskB] = maskAC;
    	bitMap[maskC] = maskAB;
    	bitMap[maskA|maskB] = maskC;
    	bitMap[maskA|maskC] = maskB;
    	bitMap[maskB|maskC] = maskA;
    	bitMap[maskA|maskB|maskC] = 0;
        for (Node node1 : g.nodes()) {
        	int id1 = node1.getId();
        	for (Node node2 : g.nodes()) {
        		int id2 = node2.getId();
        		if (node1.getId() >= node2.getId())continue;
        		if (!g.hasEdge(node1.getId(), node2.getId()))continue;
              Edge edge12 = g.getEdge(node1.getId(),node2.getId());
              if ((edge12.getColor() == capFACBond)||(edge12.getColor() == mCapFACBond)){
            	  bits[id1] = (bits[id1]  | maskA);//bits1=bits+maskA(4)
              } else if ((edge12.getColor() == capFBCBond)||(edge12.getColor() == mCapFBCBond)){
            	  bits[id1]  = (bits[id1]  | maskB);//bits1=bits+maskB(2)

              }	else if (edge12.getColor() == capFCABond||edge12.getColor() == capFCBBond||edge12.getColor() == mCapFCABond||edge12.getColor() == mCapFCBBond){
            	  bits[id1]  = (bits[id1]  | maskC);//bits1=bits+maskC(1)

              }	
              
              if (edge12.getColor() == capFCABond||edge12.getColor() == mCapFCABond){
            	  bits[id2]  = (bits[id2] | maskA);
            	  
              } else if (edge12.getColor() == capFCBBond||edge12.getColor() == mCapFCBBond){
            	  bits[id2]  = (bits[id2] | maskB);
            	  
              }	else if (edge12.getColor() == capFACBond||edge12.getColor() == capFBCBond||edge12.getColor() == mCapFACBond||edge12.getColor() == mCapFBCBond){
            	  bits[id2] = (bits[id2] | maskC);
            	  
              }		

            	}
        	char color1 = node1.getColor();
        	if (node1.getType() != Metadata.TYPE_NODE_ROOT && (foundSigmaPoint & bitMap[bits[id1]])== 0 && (color1 ==sigma0)){
        		g.getNode(node1.getId()).setColor(sigma[bits[id1]]);
        		foundSigmaPoint|= bitMap[bits[id1]];
        	}
        }
        int sigma0Factor = 0;
        int sigmaAFactor = 0;
        int sigmaBFactor = 0;
        int sigmaCFactor = 0;
        int sigmaABFactor = 0;
        int sigmaACFactor = 0;
        int sigmaBCFactor = 0;
        int sigmaABCFactor = 0;

        for(Node node1 : g.nodes()){
    	  if (node1.getType() == Metadata.TYPE_NODE_ROOT){
    		  node1.setColor('E');
    	  }
    	    else if(node1.getColor()==sigmaABCFinal){
    		    sigmaABCFactor++;
        	}
        	else if(node1.getColor()==sigmaABFinal){
        		sigmaABFactor++;
        	}
        	else if(node1.getColor()==sigmaACFinal){
        		sigmaACFactor++;
        	}
        	else if(node1.getColor()==sigmaBCFinal){
        		sigmaBCFactor++;
        	}
        	else if(node1.getColor()==sigmaAFinal){
        		sigmaAFactor++;
        	}
        	else if(node1.getColor()==sigmaBFinal){
        		sigmaBFactor++;
        	}
        	else if(node1.getColor()==sigmaCFinal){
        		sigmaCFactor++;
        	}
        	else if(node1.getColor()==sigma0Final||node1.getColor()==sigma0){
        		sigma0Factor++;
        	}
        }

          if(foundSigmaPoint == 0){
        	  if(g.factors().length != 8) {
        		  g.setNumFactors(8);
        		  g.addFactors(new int[]{g.nodeCount(),0,0,0,0,0,0,0});
        	  }
              result.add(g);
              continue;
          }

          g.setNumFactors(8);
          g.addFactors(new int[]{sigma0Factor,sigmaAFactor,sigmaBFactor,sigmaCFactor,sigmaABFactor,sigmaACFactor,sigmaBCFactor,sigmaABCFactor});//decorate only 1 point at one time
          product.clear();
          product.add(g);
          if ((foundSigmaPoint & maskA) != 0){
	          product = mulFlex.apply(product, rho0CASets, mfp);
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
              Set<Graph> product1 = new GraphList(null);
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CBSets, mfp);
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
          if ((foundSigmaPoint & maskC) != 0){
        	  Set<Graph> g1 = new HashSet<Graph>();//single diagram
              Set<Graph> product1 = new GraphList(null);
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CCSets, mfp);
    	          for (Graph gP2:newProduct){
    	              for(Node node1 : gP2.nodes()){
    	                  if(node1.getColor() == sigmaC){
    	                      if(gP2.nodeCount()== gP.nodeCount())
    	                          node1.setColor(sigmaCFinal);//sigmaC
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
              Set<Graph> product1 = new GraphList(null);
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CABSets, mfp);
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
          if ((foundSigmaPoint & maskAC) != 0){
        	  Set<Graph> g1 = new HashSet<Graph>();//single diagram
              Set<Graph> product1 = new GraphList(null);
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CACSets, mfp);
    	          for (Graph gP2:newProduct){
    	              for(Node node1 : gP2.nodes()){
    	                  if(node1.getColor() == sigmaAC){
    	                      if(gP2.nodeCount()== gP.nodeCount())
    	                          node1.setColor(sigmaACFinal);//sigmaAC
    	                      else node1.setColor(sigma0);//sigma0
    	                  }
    	              }
    	          }
    	          product1.addAll(newProduct);
              }
              product = product1;
          }
          if ((foundSigmaPoint & maskBC) != 0){
        	  Set<Graph> g1 = new HashSet<Graph>();//single diagram
              Set<Graph> product1 = new GraphList(null);
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CBCSets, mfp);
    	          for (Graph gP2:newProduct){
    	              for(Node node1 : gP2.nodes()){
    	                  if(node1.getColor() == sigmaBC){
    	                      if(gP2.nodeCount()== gP.nodeCount())
    	                          node1.setColor(sigmaBCFinal);//sigmaBC
    	                      else node1.setColor(sigma0);//sigma0
    	                  }
    	              }
    	          }
    	          product1.addAll(newProduct);
              }
              product = product1;
          }
          if ((foundSigmaPoint & maskABC) != 0){
        	  Set<Graph> g1 = new HashSet<Graph>();//single diagram
              Set<Graph> product1 = new GraphList(null);
              for (Graph gP:product){
            	  g1.clear();
            	  g1.add(gP);
            	  Set<Graph> newProduct =  mulFlex.apply(g1, rho0CABCSets, mfp);
    	          for (Graph gP2:newProduct){
    	              for(Node node1 : gP2.nodes()){
    	                  if(node1.getColor() == sigmaABC){
    	                      if(gP2.nodeCount()== gP.nodeCount())
    	                          node1.setColor(sigmaABCFinal);//sigmaABC
    	                      else node1.setColor(sigma0);//sigma0
    	                  }
    	              }
    	          }
    	          product1.addAll(newProduct);
              }
              product = product1;
          }
          pool.addAll(product);
          if (false) {
            for (Graph gP:product){
          	  int sigma0Count = checkGraph(gP);
          	  if(sigma0Count > -1){
          	    if (g.nodeCount() == gP.nodeCount()&&sigma0Count == 1){
          	      // System.out.println("@#$");
          	      coefficient.add(gP.coefficient());
          	    }
          	    //System.out.println("g "+g);
          	    //System.out.println("gP "+gP);
          	    //break;
          	  }
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
  public static void checkSet(Set<Graph> setG){
	  for (Graph g:setG){
		  if(checkGraph(g) > -1){
			  //System.out.println(g);
		  }
	  }
  }
  public static int checkGraph(Graph gP){
	  if (gP.nodeCount() != gP.edgeCount()+1 || DecorateWertheim2Site.graphHasEdgeColor(gP, 'B')){
		  return -1;
	  }
	  if (gP.nodeCount() !=4){
		  return -1;
	  }

			  boolean branch = false;
			  int sigma0Count = 0;
		  for (Node node: gP.nodes()){
			  if (node.getColor()=='A')sigma0Count++;
			  if (gP.getOutDegree(node.getId())==3){//getOutDegree: # of bonds for one point
				  branch = true;
				  byte id0 = gP.getOutNode(node.getId(), (byte)0);
				  byte id1 = gP.getOutNode(node.getId(), (byte)1);
				  byte id2 = gP.getOutNode(node.getId(), (byte)2);
				  //boolean sameDirection = (id0-node.getId())*(id1-node.getId())>0;//where the given bond is going to.. index = number of bond
				  char color0 = gP.getEdge(node.getId(), id0).getColor();
				  char color1 = gP.getEdge(node.getId(), id1).getColor();
				  char color2 = gP.getEdge(node.getId(), id2).getColor();
				  char site0 = color0=='F'?'A':(color0 == 'I'?'B':'C');//if color0 = F ->site0 is A, if not, (if color0 is I, site0 is B.. if not, site0 is C)
				  char site1 = color1=='F'?'A':(color1 == 'I'?'B':'C');
				  char site2 = color2=='F'?'A':(color2 == 'I'?'B':'C');
				  boolean sameColor = (site0 == site1)||(site0 == site2)||(site1 == site2);
				  if(sameColor){
					  return -1;
				  }
			  }
		  }
		  return branch?sigma0Count:-1;
          //break;
  }
  
  protected Set<Graph>[] splitSets(Set<Graph> graphSet, int n) {
    Set<Graph>[] sets = new Set[n+1];
    for (int i=0; i<n+1; i++) {
      sets[i] = new HashSet<Graph>();
    }
    for (Graph g : graphSet) {
      int numField2 = NumFieldNodes.value(g);
      if (numField2 < sets.length) {
        sets[numField2].add(g);
      }
    }
    return sets;
  }

  public static class DecorateWertheimParameters3Site implements Parameters {
    public final MulFlexibleParameters mfp;
    public final char capFACBond, capFBCBond, capFCABond,capFCBBond, mCapFACBond, mCapFBCBond, mCapFCABond,mCapFCBBond;
    public final Set<Graph> rhoA,rhoB,rhoC,rhoAB,rhoAC,rhoBC,rhoABC;
    public DecorateWertheimParameters3Site(MulFlexibleParameters mfp, char capFACBond, char capFBCBond, char capFCABond, char capFCBBond,char mCapFACBond, char mCapFBCBond, char mCapFCABond, char mCapFCBBond, Set<Graph> rhoA, Set<Graph> rhoB, Set<Graph> rhoC, Set<Graph> rhoAB, Set<Graph> rhoAC, Set<Graph> rhoBC, Set<Graph> rhoABC) {
      this.capFACBond = capFACBond;
      this.capFBCBond = capFBCBond;
      this.capFCABond = capFCABond;
      this.capFCBBond = capFCBBond; 
      this.mCapFACBond = mCapFACBond;
      this.mCapFBCBond = mCapFBCBond;
      this.mCapFCABond = mCapFCABond;
      this.mCapFCBBond = mCapFCBBond; 
      this.rhoA = rhoA;
      this.rhoB = rhoB;
      this.rhoC = rhoC;
      this.rhoAB = rhoAB;
      this.rhoAC = rhoAC;
      this.rhoBC = rhoBC;
      this.rhoABC = rhoABC;
      this.mfp = mfp;
    }
   
  }
}