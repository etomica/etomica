package etomica.modules.chainequilibrium;

import java.io.Serializable;

import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IBox;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.data.Data;
import etomica.data.DataSource;
import etomica.data.DataSourceIndependent;
import etomica.data.DataTag;
import etomica.data.IDataInfo;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.Null;
import etomica.units.Quantity;

/**
 * @author Matt Moynihan MoleuclarCount returns an array with the number of
 *         atoms In molecules with [1,2,3,4,5,6,7-10,10-13,13-25, <25] atoms
 */
public class MeterChainLength implements DataSource, Serializable, AgentSource, DataSourceIndependent {

    public MeterChainLength(AtomLeafAgentManager aam) {
        tag = new DataTag();
        setupData(1);
        agentManager = aam;
    }
    
    public DataTag getTag() {
        return tag;
    }

    /**
     * Creates the data object (a DataFunction) to be returned by getData().  
     * data wraps the histogram's double[] so copying is not needed.
     */
    protected void setupData(int maxChainLength) {

        xData = new DataDoubleArray(maxChainLength);
        xDataInfo = new DataInfoDoubleArray("Chain Length", Quantity.DIMENSION, new int[]{maxChainLength});
        xDataInfo.addTag(tag);
        double[] x = xData.getData();
        for (int i=0; i<maxChainLength; i++) {
            x[i] = i+1;
        }

        DataFunction newData = new DataFunction(new int[]{maxChainLength});
        if (data != null) {
            double[] y = newData.getData();
            for (int i=0; i<data.getLength() && i<y.length; i++) {
                y[i] = data.getValue(i);
            }
        }
        data = newData;
        dataInfo = new DataInfoFunction("Chain Length Distribution", Null.DIMENSION, this);
        dataInfo.addTag(tag);
    }
    
    public Class getAgentClass() {
        return AtomTag.class;
    }
    
    public Object makeAgent(IAtomLeaf a) {
        return new AtomTag();
    }
    
    // does nothing
    public void releaseAgent(Object agent, IAtomLeaf atom) {}

    //returns the number of molecules with [1,2,3,4,5,6,7-10,10-13,13-25, >25]
    // atoms
    public Data getData() {
        
        double[] histogram = data.getData();
        for (int i=0; i<histogram.length; i++) {
            histogram[i] = 0;
        }

        // untag all the Atoms
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            ((AtomTag)tagManager.getAgent(leafList.getAtom(i))).tagged = false;
        }

        for (int i=0; i<nLeaf; i++) {
            IAtomLeaf a = (IAtomLeaf)leafList.getAtom(i);
            // if an Atom is tagged, it was already counted as part of 
            // another chain
            if (((AtomTag)tagManager.getAgent(a)).tagged) continue;

            int chainLength = recursiveTag(a);
            
            if (chainLength > histogram.length) {
                setupData(chainLength);
                histogram = data.getData();
            }
            histogram[chainLength-1] += chainLength;
        }

        for (int i=0; i<histogram.length; i++) {
            histogram[i] /= box.getLeafList().getAtomCount();
        }
        
        return data;
    }

    public void reset() {
        setupData(1);
    }
    
    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }
    
    public int getIndependentArrayDimension() {
        return 1;
    }

    protected int recursiveTag(IAtomLeaf a) {
        ((AtomTag)tagManager.getAgent(a)).tagged = true;

        IAtomLeaf[] nbrs = (IAtomLeaf[])agentManager.getAgent(a);

        int ctr = 1;
        
        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(((AtomTag)tagManager.getAgent(nbrs[i])).tagged) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            ctr += recursiveTag(nbrs[i]);
        }
        return ctr;
        
    }
    
    public IBox getBox() {
        return box;
    }
    
    public void setBox(IBox box) {
        this.box = box;
        if (tagManager != null) {
            // allow old agentManager to de-register itself as a BoxListener
            tagManager.dispose();
        }
        tagManager = new AtomLeafAgentManager(this,box);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private AtomLeafAgentManager tagManager;
    private AtomLeafAgentManager agentManager;
    private DataFunction data;
    private DataDoubleArray xData;
    private DataInfoDoubleArray xDataInfo;
    private DataInfoFunction dataInfo;
    private final DataTag tag;
    
    public static class AtomTag {
        public boolean tagged;
    }

}
