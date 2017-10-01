/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataSource;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

import java.util.Arrays;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirialExternalFieldSW implements IDataSource, java.io.Serializable {

    private double lambdaWF;
	private double temperature;
	private double epsilonWF;

	/**
	 * Constructor for MeterVirial.
	 */
	public MeterVirialExternalFieldSW(ClusterAbstract aCluster, double[] wallposition, double lambdaWF, double temperature, double epsilonWF) {
		cluster = aCluster;
		this.wallPosition = wallposition;
		this.lambdaWF = lambdaWF;
		this.temperature =temperature;
		this.epsilonWF =epsilonWF;
        data = new DataDoubleArray(wallposition.length+1);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{wallposition.length+1});
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
    	IAtomList atoms = box.getLeafList();
    	double [] z = new double [atoms.getAtomCount()+1];
    	for (int i=0; i<atoms.getAtomCount()+1;i++){
    		if (i == atoms.getAtomCount()){
    			z[i]=z[i-1]+2;
    			break;
    		}
        	z[i] = atoms.getAtom(i).getPosition().getX(2);
        	
    	}
    	Arrays.sort(z);
    	
        double pi = box.getSampleCluster().value(box);
        double x[] = data.getData();
        double v=cluster.value(box)/pi;
        
   
        for (int i=0; i<wallPosition.length; i++) {
            if (z[0]-0.5 < wallPosition[i]){
            	x[i] = 0;
            }
            else x[i]=v;
        }
		double a = wallPosition[0];
		x[x.length - 1] = 0;
		double c = z[0] - 0.5;
		for (int i = 0; i < atoms.getAtomCount()+1; i++) {
			double b = z[i] - (1 + lambdaWF) * 0.5;
			if (b > c) {
				x[x.length - 1] += (c - a) * Math.exp(i*epsilonWF / temperature) * v;
				break;
			}
			x[x.length - 1] += (b - a) * Math.exp(i*epsilonWF / temperature) * v;
			a = b;
		}
    	  
        return data;
    }
    
    public ClusterAbstract getCluster() {
        return cluster;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final ClusterAbstract cluster;
	private final DataDoubleArray data;
	private final IEtomicaDataInfo dataInfo;
    private final DataTag tag;
    private final double [] wallPosition;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
    
    protected static class Near implements Comparable<Near>{
    	public IAtom atom;

		public int compareTo(Near o) {
			if (atom.getPosition().getX(2) > o.atom.getPosition().getX(2)){
				return 1;
			}else if (atom.getPosition().getX(2) < o.atom.getPosition().getX(2)){
				return -1;
			}
				else{
					return 0;
				}
			}
		}
}

