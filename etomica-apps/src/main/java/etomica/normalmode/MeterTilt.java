/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
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

    public MeterTilt(Space space, ISpecies species, int nPlanes) {
        this.species = species;
        dr = space.makeVector();
        drSum = new Vector[nPlanes+1];
        for (int i=0; i<nPlanes+1; i++) {
            drSum[i] = space.makeVector();
        }
        data = new DataDoubleArray(nPlanes+1);
        dataInfo = new DataInfoDoubleArray("tilt", Angle.DIMENSION, new int[]{nPlanes+1});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }

    public IData getData() {
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.size();
        for (int i=0; i<drSum.length; i++) {
            drSum[i].E(0);
        }
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.get(i);
            IAtomList atomList = molecule.getChildList();
            int leafCount = atomList.size();
            dr.E(atomList.get(leafCount-1).getPosition());
            dr.ME(atomList.get(0).getPosition());
            drSum[0].PE(dr);
            int iPlane = (i/2)%(drSum.length-1);
            drSum[iPlane+1].PE(dr);
        }
        double[] x = data.getData();
        for (int i=0; i<x.length; i++) {
            x[i] = Math.acos(drSum[i].getX(2)/Math.sqrt(drSum[i].squared()));
        }
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
    protected final Vector dr;
    protected final Vector[] drSum;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
}
