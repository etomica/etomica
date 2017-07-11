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
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.graph.model.Graph;
import etomica.graph.property.HasSimpleArticulationPoint;
import etomica.graph.traversal.BCVisitor;
import etomica.units.dimensions.Null;
import etomica.virial.cluster.ExternalVirialDiagrams;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirialExternalFieldRho implements IEtomicaDataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialExternalFieldRho(ExternalVirialDiagrams diagrams, MayerFunction f, double[] wallposition) {
		Set<Graph> gset = diagrams.getMSMCGraphsEX(true);
		clusters = new ArrayList<ClusterAbstract>(gset.size());
		listsPoint = new ArrayList<List<List<Byte>>>(gset.size());
		
		/*int gNum = 0;*/
		for(Graph g : gset){
			/*if (gNum < 2) {
				gNum++;continue;
			}
			gNum++;*/
			ArrayList<ClusterBonds> allBonds = new ArrayList<ClusterBonds>();
            ArrayList<Double> weights = new ArrayList<Double>();
            diagrams.populateEFBonds(g, allBonds, weights, false);  
            double [] w = new double[]{weights.get(0)};            
            
            clusters.add(new ClusterSum(allBonds.toArray(new ClusterBonds[0]), w, new MayerFunction[]{f}));
            ArrayList<List<Byte>> listComponent = new ArrayList<List<Byte>>();
            ArrayList<Byte> listPointG = new ArrayList<Byte>();
            listComponent.add(listPointG);
            List<List<Byte>> biComponents = BCVisitor.getBiComponents(g);
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
            /*System.out.println("g="+g);
            break;*/
		}
		for (ClusterAbstract cluster : clusters){
			cluster.setTemperature(1);
		}
		this.wallPosition = wallposition;
        data = new DataDoubleArray(wallposition.length+2);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{wallposition.length+2});
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
    	double pi = box.getSampleCluster().value(box);
        double x[] = data.getData();
        IAtomList atoms = box.getLeafList();
        /*double lowestatom = 0;
        for (int i=1; i<atoms.getAtomCount();i++){
        	double z = atoms.getAtom(i).getPosition().getX(2);
        	if (z<lowestatom){
        		lowestatom = z;      		
        	}
        }*/
        
        x[0]=0;
        
        for (int i=0; i<wallPosition.length; i++) {
        	double vSum = 0;
        	for(int c=0; c<clusters.size();c++){
        		double v =clusters.get(c).value(box)/pi;
        		 if (v==0){
         	    	continue;
         	    }
        		if(i==0) x[0]+=v;
        	    List<List<Byte>> cList = listsPoint.get(c);
        	    List<Byte> gList = cList.get(0);
        	    for(byte g : gList){
        	    	if (atoms.getAtom(g).getPosition().getX(2)-0.5 < wallPosition[i]){
        	    		v = 0;
        	    		break;
        	    	}
        	    }
        	    if (v==0){
        	    	continue;
        	    }
        	    for(int j=1; j<cList.size(); j++){
        	    	List<Byte> gm1List = cList.get(j);
        	    	boolean inwall = false;
        	    	for(byte gm1 :gm1List){
        	    		if (atoms.getAtom(gm1).getPosition().getX(2)-0.5 < wallPosition[i]){
        	    			inwall = true;
        	    			break;
        	    		}
        	    	}
        	    	if(!inwall){
        	    		v = 0;
        	    		break;
        	    	
        	    	}
        	    	v=-v;
        	    }
        	    vSum += v;
            }
        	x[i+1]=vSum;
        }
        /*x[x.length-1]=(atoms.getAtomCount()-1+lowestatom)*v;*/
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
