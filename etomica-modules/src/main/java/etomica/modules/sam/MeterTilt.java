/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDouble;
import etomica.data.types.DataGroup;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Angle;

/**
 * Meter that measures the average tilt angle (not the angle of average tilt!)
 *
 * @author Andrew Schultz
 */
public class MeterTilt implements IDataSource {

    public MeterTilt(Space space, ISpecies species) {
        this.species = species;
        dr = space.makeVector();
        drSum = space.makeVector();
        tag = new DataTag();
        data = new DataGroup(new IData[]{new DataDouble(), new DataDouble()});
        dataInfo = new DataGroup.DataInfoGroup("Tilt", Angle.DIMENSION, new IDataInfo[]{
                new DataDouble.DataInfoDouble("Tilt", Angle.DIMENSION),
                new DataDouble.DataInfoDouble("Tilt", Angle.DIMENSION)});
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }

    public IData getData() {
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.getMoleculeCount();
        drSum.E(0);
        double thetaSum = 0;
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atomList = molecule.getChildList();
            int leafCount = atomList.size();
            dr.E(atomList.get(leafCount-1).getPosition());
            dr.ME(atomList.get(1).getPosition());
            thetaSum += Math.acos(dr.getX(1)/Math.sqrt(dr.squared()));
            drSum.PE(dr);
        }
        ((DataDouble)data.getData(0)).x = Math.acos(drSum.getX(1)/Math.sqrt(drSum.squared()));
        ((DataDouble)data.getData(1)).x = thetaSum / nMolecules;
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    private static final long serialVersionUID = 1L;
    protected final ISpecies species;
    protected Box box;
    protected final Vector dr, drSum;
    protected final DataTag tag;
    protected final DataGroup data;
    protected final IDataInfo dataInfo;
}
