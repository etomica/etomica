package etomica.virial;

import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.units.Null;

/**
 * Measures value of clusters in a box and returns the values
 * divided by the sampling bias from the sampling cluster.
 */
public class MeterVirial implements DataSource, java.io.Serializable {

    /**
	 * Constructor for MeterVirial.
	 */
	public MeterVirial(ClusterAbstract[] aClusters) {
		clusters = aClusters;
        data = new DataDoubleArray(clusters.length);
        dataInfo = new DataInfoDoubleArray("Cluster Value",Null.DIMENSION, new int[]{clusters.length});
        tag = new DataTag();
        dataInfo.addTag(tag);
	}

	public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public Data getData() {
        CoordinatePairSet cPairSet = box.getCPairSet();
        AtomPairSet aPairSet = box.getAPairSet();
        double pi = box.getSampleCluster().value(cPairSet,aPairSet);
        double x[] = data.getData();
        for (int i=0; i<clusters.length; i++) {
        	
    		if (clusters[i] instanceof ClusterCoupledFlipped) {
    			((ClusterCoupledFlipped)clusters[i]).setPhase(box);
    		}
    		else if (clusters[i] instanceof ClusterWeightAbs) {
    		    ClusterAbstract weightCluster = ((ClusterWeightAbs)clusters[i]).getSubCluster();
    		    if (weightCluster instanceof ClusterCoupledFlipped) {
    		        ((ClusterCoupledFlipped)weightCluster).setPhase(box);
    		    }
    		}
        	
            x[i] = clusters[i].value(cPairSet,aPairSet)/pi;
            if (Double.isNaN(x[i])) throw new RuntimeException();
        }
        return data;
    }
    
    public BoxCluster getBox() {
        return box;
    }
    
    public void setBox(BoxCluster newBox) {
        box = newBox;
    }

    protected final ClusterAbstract clusters[];
	private final DataDoubleArray data;
	private final IDataInfo dataInfo;
    private final DataTag tag;
    private BoxCluster box;
    private static final long serialVersionUID = 1L;
}
