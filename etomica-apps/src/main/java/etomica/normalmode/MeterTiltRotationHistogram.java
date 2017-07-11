/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramSimple;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.data.types.DataGroup;
import etomica.data.types.DataGroup.DataInfoGroup;
import etomica.math.DoubleRange;
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
public class MeterTiltRotationHistogram implements IAction, IEtomicaDataSource {

    public MeterTiltRotationHistogram(Space space, ISpecies species, int nPlanes) {
        this.species = species;
        int nData = 360;
        dr = space.makeVector();
        drSum = new Vector[nPlanes];
        for (int i=0; i<nPlanes; i++) {
            drSum[i] = space.makeVector();
        }
        histogram = new HistogramSimple[2];
        double[] lim = new double[]{Math.PI, Math.PI/4};
        xDataSource = new DataSourceIndependentSimple[2];
        data = new DataFunction[2];
        iTag = new DataTag[2];
        dataInfo = new DataInfoFunction[2];
        for (int i=0; i<2; i++) {
            histogram[i] = new HistogramSimple(nData, new DoubleRange(-lim[i], lim[i]));
            DataInfoDoubleArray xDataInfoi = new DataInfoDoubleArray("angle", Angle.DIMENSION, new int[]{nData});
            DataTag xTagi = new DataTag();
            xDataInfoi.addTag(xTagi);
            xDataSource[i] = new DataSourceIndependentSimple(histogram[i].xValues(), xDataInfoi);
            data[i] = new DataFunction(new int[]{nData}, histogram[i].getHistogram());
            dataInfo[i] = new DataInfoFunction("tilt", Angle.DIMENSION, xDataSource[i]);
            iTag[i] = new DataTag();
            dataInfo[i].addTag(iTag[i]);
        }
        dataGroup = new DataGroup(data);
        dataInfoGroup = new DataInfoGroup("tilt histograms", Angle.DIMENSION, dataInfo);
        groupTag = new DataTag();
        dataInfoGroup.addTag(groupTag);
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public void reset() {
        for (int i=0; i<histogram.length; i++) {
            histogram[i].reset();
        }
    }

    public void actionPerformed() {
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.getMoleculeCount();
        int nPlanes = drSum.length;
        for (int i=0; i<nPlanes; i++) {
            drSum[i].E(0);
        }
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atomList = molecule.getChildList();
            int leafCount = atomList.getAtomCount();
            dr.E(atomList.getAtom(leafCount-1).getPosition());
            dr.ME(atomList.getAtom(0).getPosition());
            histogram[0].addValue(Math.atan2(dr.getX(1), dr.getX(0)));
            int iPlane = (i/2)%nPlanes;
            drSum[iPlane].PE(dr);
        }
        for (int i=0; i<nPlanes; i++) {
            histogram[1].addValue(Math.atan2(drSum[i].getX(1), drSum[i].getX(0)));
        }
    }

    public IData getData() {
        for (int i=0; i<histogram.length; i++) {
            histogram[i].getHistogram();
        }
        return dataGroup;
    }

    public IEtomicaDataInfo getDataInfo() {
        return dataInfoGroup;
    }

    public DataTag getTag() {
        return groupTag;
    }
    
    public DataInfoFunction getIDataInfo(int i) {
        return dataInfo[i];
    }
    
    public DataTag getTag(int i) {
        return iTag[i];
    }

    private static final long serialVersionUID = 1L;
    protected final ISpecies species;
    protected Box box;
    protected final Vector dr;
    protected final Vector[] drSum;
    protected final DataGroup dataGroup;
    protected final DataInfoGroup dataInfoGroup;
    protected final DataFunction[] data;
    protected final DataInfoFunction[] dataInfo;
    protected final DataSourceIndependentSimple[] xDataSource;
    protected final DataTag groupTag;
    protected final DataTag[] iTag;
    protected final HistogramSimple[] histogram;
}
