/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISpecies;
import etomica.api.IVectorMutable;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IEtomicaDataInfo;
import etomica.data.IEtomicaDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.space.ISpace;
import etomica.units.Angle;

/**
 * Meter that measures the average tilt angle (not the angle of average tilt!)
 *
 * @author Andrew Schultz
 */
public class MeterTilt implements IEtomicaDataSource {

    public MeterTilt(ISpace space, ISpecies species, int nPlanes) {
        this.species = species;
        dr = space.makeVector();
        drSum = new IVectorMutable[nPlanes+1];
        for (int i=0; i<nPlanes+1; i++) {
            drSum[i] = space.makeVector();
        }
        data = new DataDoubleArray(nPlanes+1);
        dataInfo = new DataInfoDoubleArray("tilt", Angle.DIMENSION, new int[]{nPlanes+1});
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public void setBox(IBox newBox) {
        box = newBox;
    }

    public IData getData() {
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.getMoleculeCount();
        for (int i=0; i<drSum.length; i++) {
            drSum[i].E(0);
        }
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atomList = molecule.getChildList();
            int leafCount = atomList.getAtomCount();
            dr.E(atomList.getAtom(leafCount-1).getPosition());
            dr.ME(atomList.getAtom(0).getPosition());
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

    public IEtomicaDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    private static final long serialVersionUID = 1L;
    protected final ISpecies species;
    protected IBox box;
    protected final IVectorMutable dr;
    protected final IVectorMutable[] drSum;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
}
