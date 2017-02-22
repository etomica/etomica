/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * 
 */
package etomica.data.meter;

import etomica.api.IAtom;
import etomica.api.IAtomKinetic;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IVectorMutable;
import etomica.data.DataInfo;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataVector;
import etomica.data.types.DataVector.DataInfoVector;
import etomica.space.ISpace;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Length;
import etomica.units.Mass;
import etomica.units.Time;

/**
 * Returns the instantaneous total center-of-mass momentum, summed over all
 * leaf atoms in a box.
 *
 */
public class MeterMomentumCOM implements IEtomicaDataSource, java.io.Serializable {

    public MeterMomentumCOM(ISpace space) {
        data = new DataVector(space);
        momentumSum = data.x;
        dataInfo = new DataInfoVector("COM momentum", new CompoundDimension(
                new Dimension[] {Mass.DIMENSION, Length.DIMENSION, Time.DIMENSION}, new double[] {1.,1.,-1.}),
                space);
        tag = new DataTag();
        dataInfo.addTag(tag);
        
    }

    /**
     * Returns the instantaneous total center-of-mass momentum over all atoms in the box.
     */
    public IData getData() {
        momentumSum.E(0.0);
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtomKinetic a = (IAtomKinetic)leafList.getAtom(iLeaf);
            double mass = ((IAtom)a).getType().getMass();
            momentumSum.PEa1Tv1(mass,a.getVelocity());
        }
        return data;
    }
    
    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(IBox box) {
        this.box = box;
    }
    
    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    private static final long serialVersionUID = 1L;
    private IBox box;
    private final IVectorMutable momentumSum;
    private final DataVector data;    
    private final DataInfo dataInfo;
    private String name;
    protected final DataTag tag;

}
