/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.interfacial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceMolecular;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.types.DataDouble;
import etomica.data.types.DataDouble.DataInfoDouble;
import etomica.molecule.IMolecule;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.Angle;

/**
 * Meter for collecting the molecular orientation of the dimer.  The value
 * returned is cos(theta), where theta is the angle the dimer makes with the
 * x axis.
 */
public class MeterOrientation implements DataSourceMolecular {
    
    public MeterOrientation(Space space) {
        dataInfo = new DataInfoDouble("orientation", Angle.DIMENSION);
        data = new DataDouble();
        tag = new DataTag();
        dr = space.makeVector();
    }
    
    public DataTag getTag() {
        return tag;
    }
    
    public void setBox(Box newBox) {
        boundary = newBox.getBoundary();
    }
    
    public IData getData(IMolecule atom) {
        IAtomList children = atom.getChildList();
        dr.Ev1Mv2(children.getAtom(children.getAtomCount()-1).getPosition(),
                  children.getAtom(0).getPosition());
        boundary.nearestImage(dr);
        data.x= dr.getX(0) / Math.sqrt(dr.squared());
        return data;
    }
    
    public IEtomicaDataInfo getMoleculeDataInfo() {
        return dataInfo;
    }
    
    private static final long serialVersionUID = 1L;
    protected final DataInfoDouble dataInfo;
    protected final DataDouble data;
    protected final DataTag tag;
    protected final Vector dr;
    protected Boundary boundary;
}
