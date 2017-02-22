/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/**
 * 
 */
package etomica.data.meter;

import etomica.api.IAtom;
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
import etomica.units.Length;

/**
 * Returns the instantaneous center-of-mass position, summed over all
 * leaf atoms in a box, dividing by the number of atoms.
 *
 */
public class MeterPositionCOM implements IEtomicaDataSource, java.io.Serializable {

    public MeterPositionCOM(ISpace space) {
        data = new DataVector(space);
        positionSum = data.x;
        dataInfo = new DataInfoVector("COM momentum", Length.DIMENSION, space);
        tag = new DataTag();
        dataInfo.addTag(tag);
        
    }

    /**
     * Returns the position of the center of mass of all atoms in the box.
     */
    public IData getData() {
        positionSum.E(0.0);
        double massSum = 0.0;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            IAtom a = leafList.getAtom(iLeaf);
            double mass = a.getType().getMass();
            massSum += mass;
            positionSum.PEa1Tv1(mass,a.getPosition());
        }
        positionSum.TE(1.0/massSum);
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
    public void setBox(IBox newBox) {
        box = newBox;
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
    private final IVectorMutable positionSum;
    private final DataVector data;    
    private final DataInfo dataInfo;
    private String name;
    protected final DataTag tag;

}
