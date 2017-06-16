/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.data.DataSourceAtomic;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataTensor;
import etomica.data.types.DataTensor.DataInfoTensor;
import etomica.space.Space;
import etomica.units.Energy;

/**
 * A meter to compute the velocity component of the pressure tensor. 
 * Averages a tensor quantity formed from a dyad of the velocity of each atom. 
 * Specifically, the quantity averaged is 1/N * sum(pp/m), where p is the momentum,
 * m is the mass, and the sum is over all N atoms.
 * 
 * @author Rob Riggleman
 */

public class MeterTensorVelocity implements DataSourceAtomic, java.io.Serializable {

    public MeterTensorVelocity(Space space) {
        data = new DataTensor(space);
        dataInfo = new DataInfoTensor("pp/m",Energy.DIMENSION, space);
        atomData = new DataTensor(space);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }
    
    public DataTag getTag() {
        return tag;
    }
       
    public IEtomicaDataInfo getAtomDataInfo() {
        return dataInfo;
    }
       
    /**
     * Returns the velocity dyad (mass*vv) summed over all atoms, and divided by N
     */
    public IData getData() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        data.E(0.0);
        int count = 0;
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
            getData(leafList.getAtom(iLeaf));
            data.PE(atomData);
            count++;
        }
        data.TE(1.0/nLeaf);
        return data;
    }
    
    /**
     * Returns the velocity dyad (mass*vv) for the given atom.
     */
    public IData getData(IAtom atom) {
        Vector vel = ((IAtomKinetic)atom).getVelocity();
        atomData.x.Ev1v2(vel, vel);
        atomData.TE(atom.getType().rm());
        return atomData;
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }
    /**
     * @param box The box to set.
     */
    public void setBox(Box box) {
        this.box = box;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }
    
    private static final long serialVersionUID = 1L;
    private String name;
    private Box box;
    private final DataTensor data, atomData;
    private final DataInfoTensor dataInfo;
    protected DataTag tag;
}
