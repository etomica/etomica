/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.units.dimensions.Null;
import etomica.units.dimensions.Quantity;

import java.io.Serializable;

/**
 * @author Matt Moynihan MoleuclarCount returns an array with the number of
 *         atoms In molecules with [1,2,3,4,5,6,7-10,10-13,13-25, <25] atoms
 */
public class MeterChainLength implements IDataSource, Serializable, AgentSource<MeterChainLength.AtomTag>, DataSourceIndependent {

    private static final long serialVersionUID = 1L;
    protected final DataTag tag, xTag;
    protected Box box;
    protected AtomLeafAgentManager<AtomTag> tagManager;
    protected AtomLeafAgentManager<IAtom[]> agentManager;
    protected DataFunction data;
    protected DataDoubleArray xData;
    protected DataInfoDoubleArray xDataInfo;
    protected DataInfoFunction dataInfo;
    protected AtomType ignoredAtomType;

    public MeterChainLength(AtomLeafAgentManager<IAtom[]> aam) {
        tag = new DataTag();
        xTag = new DataTag();
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
        xDataInfo.addTag(xTag);
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

    public AtomTag makeAgent(IAtom a, Box agentBox) {
        return new AtomTag();
    }

    // does nothing
    public void releaseAgent(AtomTag agent, IAtom atom, Box agentBox) {}

    //returns the number of molecules with [1,2,3,4,5,6,7-10,10-13,13-25, >25]
    // atoms
    public IData getData() {

        double[] histogram = data.getData();
        for (int i=0; i<histogram.length; i++) {
            histogram[i] = 0;
        }

        // untag all the Atoms
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            tagManager.getAgent(leafList.getAtom(i)).tagged = false;
        }

        int totalAtoms = 0;
        for (int i=0; i<nLeaf; i++) {
            IAtom a = leafList.getAtom(i);
            if (a.getType() == ignoredAtomType) continue;
            // if an Atom is tagged, it was already counted as part of
            // another chain
            if (tagManager.getAgent(a).tagged) continue;

            int chainLength = recursiveTag(a);

            if (chainLength > histogram.length) {
                setupData(chainLength);
                histogram = data.getData();
            }
            histogram[chainLength-1] += chainLength;
            totalAtoms += chainLength;
        }

        for (int i=0; i<histogram.length; i++) {
            histogram[i] /= totalAtoms;
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

    public DataTag getIndependentTag() {
        return xTag;
    }

    protected int recursiveTag(IAtom a) {
        if (a.getType() == ignoredAtomType) return 0;
        tagManager.getAgent(a).tagged = true;

        IAtom[] nbrs = agentManager.getAgent(a);

        int ctr = 1;

        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if (tagManager.getAgent(nbrs[i]).tagged) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            ctr += recursiveTag(nbrs[i]);
        }
        return ctr;

    }

    public Box getBox() {
        return box;
    }

    public void setBox(Box box) {
        this.box = box;
        if (tagManager != null) {
            // allow old agentManager to de-register itself as a BoxListener
            tagManager.dispose();
        }
        tagManager = new AtomLeafAgentManager<AtomTag>(this,box);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public AtomType getIgnoredAtomType() {
        return ignoredAtomType;
    }

    public void setIgnoredAtomType(AtomType ignoredAtomType) {
        this.ignoredAtomType = ignoredAtomType;
    }
    
    public static class AtomTag {
        public boolean tagged;
    }

}
