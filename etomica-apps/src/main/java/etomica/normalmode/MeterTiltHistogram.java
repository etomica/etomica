/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.action.IAction;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.histogram.HistogramNotSoSimple;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
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
public class MeterTiltHistogram implements IAction, IDataSource, DataSourceIndependent {

    public MeterTiltHistogram(Space space, ISpecies species) {
        this.species = species;
        int nData = 180;
        dr = space.makeVector();
        histogram = new HistogramNotSoSimple(nData, new DoubleRange(0, Math.PI/2));
        histogram.setDoAveraging(false);
        xData = new DataDoubleArray(new int[]{nData}, histogram.xValues());
        xDataInfo = new DataInfoDoubleArray("angle", Angle.DIMENSION, new int[]{nData});
        data = new DataFunction(new int[]{nData}, histogram.getHistogram());
        dataInfo = new DataInfoFunction("tilt", Angle.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
        xTag = new DataTag();
        xDataInfo.addTag(xTag);
    }
    
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    public void reset() {
        histogram.reset();
    }

    public void actionPerformed() {
        IMoleculeList molecules = box.getMoleculeList(species);
        int nMolecules = molecules.getMoleculeCount();
        for (int i=0; i<nMolecules; i++) {
            IMolecule molecule = molecules.getMolecule(i);
            IAtomList atomList = molecule.getChildList();
            int leafCount = atomList.getAtomCount();
            dr.E(atomList.getAtom(leafCount-1).getPosition());
            dr.ME(atomList.getAtom(0).getPosition());
            dr.normalize();
            double sintheta = Math.sqrt(dr.getX(0)*dr.getX(0) + dr.getX(1)*dr.getX(1));
            double costheta = dr.getX(2);
            double theta = Math.atan2(sintheta, costheta);
            histogram.addValue(theta, 1/sintheta);
        }
    }

    public IData getData() {
        histogram.getHistogram();
        return data;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
    
    public DataDoubleArray getIndependentData(int i) {
        return xData;
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return xDataInfo;
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xTag;
    }

    private static final long serialVersionUID = 1L;
    protected final ISpecies species;
    protected Box box;
    protected final Vector dr;
    protected final DataFunction data;
    protected final DataInfoFunction dataInfo;
    protected final DataDoubleArray xData;
    protected final DataInfoDoubleArray xDataInfo;
    protected final DataTag tag, xTag;
    protected final HistogramNotSoSimple histogram;
}
