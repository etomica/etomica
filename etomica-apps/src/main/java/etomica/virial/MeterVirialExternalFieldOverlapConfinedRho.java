/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;

import etomica.atom.IAtomList;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.graph.model.Graph;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.traversal.BCVisitor;
import etomica.graph.traversal.Biconnected;
import etomica.units.Null;
import etomica.virial.cluster.ExternalVirialDiagrams;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirialExternalFieldOverlapConfinedRho implements ClusterWeightSumWall.DataSourceClusterWall, java.io.Serializable {
	
	private double walldistance;
	private double lambdaWF;
	private double temperature;
	private double epsilonWF;
    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialExternalFieldOverlapConfinedRho(ExternalVirialDiagrams diagrams, MayerFunction f, double [] wallposition, double walldistance, double lambdaWF, double temperature, double epsilonWF) {
		Set<Graph> gset = diagrams.getMSMCGraphsEX(true);
		clusters = new ArrayList<ClusterAbstract>(gset.size());
		listsPoint = new ArrayList<List<List<Byte>>>(gset.size());
	    
		//int gNum = 0;
		
		for(Graph g : gset){	
				
			//if (gNum < 5) {
			//	gNum++;continue;
			//	}	 
			 
			ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            ArrayList<Double> weights = new ArrayList<Double>();
            diagrams.populateEFBonds(g, allBonds, weights, false);  
            double [] w = new double[]{weights.get(0)};            

            clusters.add(new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f}));
            ArrayList<List<Byte>> listComponent = new ArrayList<List<Byte>>();
            ArrayList<Byte> listPointG = new ArrayList<Byte>();
            listPointG.add((byte)0);
            listComponent.add(listPointG);
            List<List<Byte>> biComponents = new ArrayList<List<Byte>>();
            BCVisitor v = new BCVisitor(biComponents);
            new Biconnected().traverseAll(g, v);
            HasSimpleArticulationPoint hap = new HasSimpleArticulationPoint();
            hap.check(g);
            
            List<List<Byte>> prvLayerList = new ArrayList<List<Byte>>();
            List<Byte> prvApList = new ArrayList<Byte>();
            prvApList.add((byte)0);
            while (!prvApList.isEmpty()){
            	List<List<Byte>> layerList = new ArrayList<List<Byte>>();
	            List<Byte> apList = new ArrayList<Byte>();
		        for(byte prvAP : prvApList){
		            
		            for(List<Byte> biComponent : biComponents){
		            	if(!biComponent.contains(prvAP) || prvLayerList.contains(biComponent)){
		            		continue;
		            	}
		        		boolean isterminal = true;
		        		
		        		for(byte b : biComponent){
		        			if(b == prvAP){
		        				continue;
		        			}
		        			if(hap.getArticulationPoints().contains(b)){
		        				isterminal = false;
		        				apList.add(b);
		        				layerList.add(biComponent);
		        					
		        			}
		        			
		        		}
		        		List<Byte> listAll = new ArrayList<Byte>();
		        		listAll.addAll(biComponent);
		        		listAll.remove((Byte)prvAP);
		        		if(isterminal){
		        			listComponent.add(listAll);
		        		}
		        		else {
		        			listPointG.addAll(listAll);
		    			}
		            }
		        }
		        prvApList = apList;
		        prvLayerList = layerList;
            }
            listsPoint.add(listComponent);
           
            //if(gNum == 5){
            //	break;
            //}
		}
		
		for (ClusterAbstract cluster : clusters){
			cluster.setTemperature(temperature);
		}
		
		this.wallPosition = wallposition; 
		this.walldistance = walldistance;
		this.lambdaWF = lambdaWF;
		this.temperature =temperature;
		this.epsilonWF =epsilonWF;
        data = new DataDoubleArray(wallposition.length + 3);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{wallposition.length + 3});
        tag = new DataTag();
        dataInfo.addTag(tag);
	}

	public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IData getData() {
    	
        double x[] = data.getData();
        IAtomList atoms = box.getLeafList();
        for(int i=0; i<x.length; i++){
        	x[i]=0;
        }          
        
        //atoms.getAtom(1).getPosition().setX(2,+0.01);
       // atoms.getAtom(2).getPosition().setX(2,+0.6);
        //atoms.getAtom(3).getPosition().setX(2,-0.6);  
            
    	for(int c=0; c<clusters.size();c++){
    		double v =clusters.get(c).value(box);       
    		//v = Math.exp(1/0.6)-1;
    		if (v==0){
    	    	continue;
    	    }
    		
    	    List<List<Byte>> cList = listsPoint.get(c);
    	    List<Byte> gList = cList.get(0);    	    
    	 
    	    double glowestatom = 0.0;
    	    double ghighestatom = 0.0;     	   	    
    	    double [] z = new double [atoms.getAtomCount()];       	   
    
    	    for(byte g : gList){
    	    	z[g] = atoms.getAtom(g).getPosition().getX(2);
    	    	if (z[g] < glowestatom){
    	    		glowestatom = z[g];    	    			
    	    	}    	 
    	    	else if (z[g] > ghighestatom){
    	    		ghighestatom = z[g];
    	    	}
    	    }  
    	   	     	  
    		x[0] += v;
    	    	
    	    	int istart = (int)Math.ceil((ghighestatom)/0.01);
    	    	int ilast = (int)Math.floor((glowestatom - 1 + walldistance)/0.01) ;
    	    	for (int i = istart; i <= ilast; i++){

    		      //for (int i=0; i <= wallPosition.length; i++){
    	    		double wallPosition = 0.01*i - walldistance + 0.5;    	    		
    	    		double g1 = 1;
    	    		for(byte g : gList){    	    			
    	    			
	    			    	//if (z[g] < wallPosition + 0.5 || wallPosition + walldistance - 0.5 < z[g]){
	    					//    g1 *= 0;
	    				    //}        	    			
	    			    	if (z[g] < wallPosition + 0.5 + lambdaWF){
    	    					g1 *= Math.exp(epsilonWF / temperature);    	    				
    	    				}
    	    				if (z[g] > wallPosition + walldistance - 0.5 - lambdaWF){
    	    					g1 *= Math.exp(epsilonWF / temperature); 
    	    				}    	    				
    	    			}    	    	
    	    		
    	    		double g4 = 1;	    		
    	    		for(int j=0; j<cList.size()-1; j++){
    	    			List<Byte> gm1List = cList.get(j+1);
    	    			double g2 = 1;    	    		
    	    			for(byte gm1 :gm1List){
    	    			
    	    				if (atoms.getAtom(gm1).getPosition().getX(2) < wallPosition + 0.5 || wallPosition + walldistance - 0.5 < atoms.getAtom(gm1).getPosition().getX(2)){
    	    					g2 *= 0;
    	    				}    	   
    	    				else if (atoms.getAtom(gm1).getPosition().getX(2) < wallPosition + 0.5 + lambdaWF){
    	    					g2 *= Math.exp(epsilonWF / temperature); 
    	    				}
    	    				if (atoms.getAtom(gm1).getPosition().getX(2) >  wallPosition + walldistance - 0.5 - lambdaWF){
    	    					g2 *= Math.exp(epsilonWF / temperature); 
    	    				}
    	    					    			
    	    			}   	    		  	    			
    	    			
    	    			double g3 = g2 - 1; 
    	    			g4 *= g3;
    	    		}   
    	    		
    	    		int numWallPosition = (int)Math.round((walldistance-1)/0.01) + 1 ;
    	    		int j = i > (numWallPosition-1)/2 ? numWallPosition - i - 1 : i;
    	    		/*if(i==127){
    	    			 System.out.println("j="+j+" numWallPosition="+numWallPosition);
    	    		}*/
    	    	
    	    		x[j + 1] += g1*g4*v;
    	    		if(j==0){
    	    			x[x.length-2] += 0.005*g1*g4*v;
    	    		}
    	    		else{
    	    		x[x.length-2] += 0.01*g1*g4*v;
    	    		}
        	    	x[x.length-1] += Math.abs(0.01*g1*g4*v);  
    	    	}    	    	
	    	    	  
    	    
      	}  
    	for(int i=0;i<Math.round((walldistance-1)/0.02);i++){
    		
    	x[i+1] = x[i+1]/2; 
    	}
        return data; 
        
    }
    
   
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }
    protected final List<List<List<Byte>>> listsPoint;
    protected final List<ClusterAbstract> clusters;
	private final DataDoubleArray data;
	private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private final double [] wallPosition;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;    
    
}
