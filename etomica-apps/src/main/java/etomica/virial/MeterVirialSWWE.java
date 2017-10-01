package etomica.virial;

import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.dimensions.Null;

import java.math.BigDecimal;

public class MeterVirialSWWE implements IDataSource {

    protected static final BigDecimal BDZERO = new BigDecimal(0);
	
	protected final int n, npairs;
	protected final ClusterWheatleyExtendSW clusters;
	private final DataDoubleArray data;
	private final IDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
	
	public MeterVirialSWWE(ClusterWheatleyExtendSW aClusters) {
		n = aClusters.pointCount();
		npairs = n *(n-1)/2;
		clusters = aClusters;
        data = new DataDoubleArray(npairs+1);
        dataInfo = new DataInfoDoubleArray("Cluster Array",Null.DIMENSION, new int[]{npairs+1});
        tag = new DataTag();
        dataInfo.addTag(tag);
	}
	
	public IDataInfo getDataInfo() {
        return dataInfo;
    }
	
	public DataTag getTag() {
		return tag;
	}
	
	public IData getData() {
		
//		IAtomList iAL =	box.getLeafList();//print coordinates info.
//    	for(int i=0; i<iAL.getAtomCount(); i++){
//    		System.out.print("   " + iAL.getAtom(i).getPosition()+",");
//    	}
//    	System.out.println();
//		CoordinatePairSet cPairs = box.getCPairSet();
//		for(int i=0; i<n; i++){//print distances between atoms.
//         	for(int j=i+1; j<n; j++){
//         		System.out.print(Math.sqrt(cPairs.getr2(i, j)));
//         	}
//    	}
//    	System.out.println();		
			
		double[] x = data.getData();//when x changes, data will change.
		  	  
		double[] v = clusters.valueArray(box);  
		double pi = box.getSampleCluster().value(box);
		if (pi == 0 || pi == Double.POSITIVE_INFINITY || Double.isNaN(pi)){ 
			throw new RuntimeException("oops "+pi);
		}
	     
		for(int i=0; i<v.length; i++){ 
			x[i] = v[i] / pi;
		}
		return data;
	}
	
	public ClusterAbstract getClusters() {
        return clusters;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }
}
