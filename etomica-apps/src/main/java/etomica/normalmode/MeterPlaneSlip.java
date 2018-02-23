/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.units.dimensions.Angle;

/**
 * Meter that measures the average tilt angle (not the angle of average tilt!)
 *
 * @author Andrew Schultz
 */
public class MeterPlaneSlip implements IDataSource {

    public MeterPlaneSlip(Space space, ISpecies species, int nPlanes, int nx, int ny) {
        this.species = species;
        pos = new MoleculePositionGeometricCenter(space);
        drSum = new Vector[nPlanes][2];
        for (int i=0; i<nPlanes; i++) {
            drSum[i][0] = space.makeVector();
            drSum[i][1] = space.makeVector();
        }
        data = new DataDoubleArray(nPlanes*2);
        dataInfo = new DataInfoDoubleArray("slip", Angle.DIMENSION, new int[]{nPlanes*2});
        tag = new DataTag();
        dataInfo.addTag(tag);
        offset0 = new double[nPlanes][2];
        this.nx = nx;
        this.ny = ny;
    }
    
    public void setBox(Box newBox) {
        box = newBox;
        
        int nPlanes = drSum.length;
        
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.size();
        for (int i=0; i<nPlanes; i++) {
            drSum[i][0].E(0);
            drSum[i][1].E(0);
        }
        for (int i=0; i<nMolecules; i++) {
            int iPlane = (i/2)%nPlanes;
            IAtomList atomList = molecules.get(i).getChildList();
            drSum[iPlane][0].PE(atomList.get(0).getPosition());
            drSum[iPlane][1].PE(atomList.get(1).getPosition());
        }
        int nMoleculesPerPlane = nMolecules/nPlanes;
        Vector ba = box.getBoundary().getEdgeVector(0);
        double a0 = ba.getX(0) / nx;
        Vector bb = box.getBoundary().getEdgeVector(1);
        double b0 = bb.getX(1) / ny;
        Vector bc = box.getBoundary().getEdgeVector(2);
        
        for (int i=0; i<nPlanes; i++) {
            drSum[i][0].TE(1.0/nMoleculesPerPlane);
            drSum[i][1].TE(1.0/nMoleculesPerPlane);
        }
        for (int i=0; i<nPlanes; i++) {
            if (i > 0) {
                offset0[i][0] = (drSum[i][0].getX(0) - drSum[i-1][1].getX(0))/a0;
                offset0[i][1] = (drSum[i][0].getX(1) - drSum[i-1][1].getX(1))/b0;
            }
            else {
                offset0[i][0] = (drSum[i][0].getX(0) - (drSum[nPlanes-1][1].getX(0) - bc.getX(0)))/a0;
                offset0[i][1] = (drSum[i][0].getX(1) - drSum[nPlanes-1][1].getX(1))/b0;
            }
        }
    }

    public IData getData() {
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.size();
        int nPlanes = drSum.length;
        for (int i=0; i<nPlanes; i++) {
            drSum[i][0].E(0);
            drSum[i][1].E(0);
        }
        for (int i=0; i<nMolecules; i++) {
            int iPlane = (i/2)%nPlanes;
            IAtomList atomList = molecules.get(i).getChildList();
            drSum[iPlane][0].PE(atomList.get(0).getPosition());
            drSum[iPlane][1].PE(atomList.get(1).getPosition());
        }
        int nMoleculesPerPlane = nMolecules/nPlanes;
        Vector ba = box.getBoundary().getEdgeVector(0);
        double a = ba.getX(0) / nx;
        Vector bb = box.getBoundary().getEdgeVector(1);
        double b = bb.getX(1) / ny;
        Vector bc = box.getBoundary().getEdgeVector(2);
        
        for (int i=0; i<nPlanes; i++) {
            drSum[i][0].TE(1.0/nMoleculesPerPlane);
            drSum[i][1].TE(1.0/nMoleculesPerPlane);
        }
        double[] x = data.getData();
        for (int i=0; i<nPlanes; i++) {
            if (i > 0) {
                x[2*i+0] = (drSum[i][0].getX(0) - drSum[i-1][1].getX(0))/a - offset0[i][0];
                x[2*i+1] = (drSum[i][0].getX(1) - drSum[i-1][1].getX(1))/b - offset0[i][1];
            }
            else {
                x[2*i+0] = (drSum[i][0].getX(0) - (drSum[nPlanes-1][1].getX(0) - bc.getX(0)))/a - offset0[i][0];
                x[2*i+1] = (drSum[i][0].getX(1) - drSum[nPlanes-1][1].getX(1))/b - offset0[i][1];
            }
        }
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected final ISpecies species;
    protected Box box;
    protected final Vector[][] drSum;
    protected final DataDoubleArray data;
    protected final DataInfoDoubleArray dataInfo;
    protected final DataTag tag;
    protected IMoleculePositionDefinition pos;
    protected final double[][] offset0;
    protected final int nx, ny;
}
