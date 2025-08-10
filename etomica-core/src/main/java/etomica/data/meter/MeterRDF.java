/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.action.IAction;
import etomica.atom.AtomPairTest;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

import java.util.Arrays;

/**
 * Meter for tabulation of the atomic radial distribution function (RDF).  The
 * meter takes data via actionPerformed and returns the average RDF via
 * getData.
 *
 * @author David Kofke
 */
public class MeterRDF implements IAction, IDataSource, DataSourceIndependent {

    protected final Space space;
    protected final DataSourceUniform xDataSource;
    protected final DataTag tag;
    protected final Vector dr;
    protected Box box;
    protected long[] gSum;
    protected DataFunction data;
    protected DataDoubleArray rData;
    protected AtomPairTest test;
    protected double xMax;
    protected long callCount;
    protected AtomType type1, type2;
    private IDataInfo dataInfo;
    private String name;
    protected boolean resetAfterData;
    protected final boolean singleSample;

	/**
	 * Creates meter with default to compute pair correlation for all
	 * leaf atoms in a box.
	 * @param space
	 */
    public MeterRDF(Space space) {
        this(space, false);
    }

    public MeterRDF(Space space, boolean singleSample) {
	    this.space = space;
        this.singleSample = singleSample;

        xDataSource = new DataSourceUniform("r", Length.DIMENSION);
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);

        rData = (DataDoubleArray)xDataSource.getData();
        data = new DataFunction(new int[] {rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);

        dr = space.makeVector();
        tag = new DataTag();
        dataInfo.addTag(tag);
    }
    
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setAtomType(AtomType type) {
        type1 = type;
        type2 = type;
    }

    public void setAtomTypes(AtomType type1, AtomType type2) {
        this.type1 = type1;
        this.type2 = type2;
    }

    /**
     * Zero's out the RDF sum tracked by this meter.
     */
    public void reset() {
        rData = (DataDoubleArray) xDataSource.getData();
        xMax = xDataSource.getXMax();
        data = new DataFunction(new int[]{rData.getLength()});
        gSum = new long[rData.getLength()];
        dataInfo = new DataInfoFunction("g(r)", Null.DIMENSION, this);
        dataInfo.addTag(tag);
        zeroData();
    }

    public void zeroData() {
        Arrays.fill(gSum, 0);
        callCount = 0;
    }

    /**
     * Takes the RDF for the current configuration of the given box.
     */
    public void actionPerformed() {
        if (rData != xDataSource.getData() ||
            data.getLength() != rData.getLength() ||
            xDataSource.getXMax() != xMax) {
            reset();
        }

        double xMaxSquared = xMax*xMax;
        Boundary boundary = box.getBoundary();
        // iterate over all pairs
        IAtomList atoms = box.getLeafList();
        for (int i=0; i<atoms.size(); i++) {
            IAtom atom1 = atoms.get(i);
            if (type1 != null && atom1.getType() != type1) continue;
            for (int j=i+1; j<atoms.size(); j++) {
                IAtom atom2 = atoms.get(j);
                if (type2 != null && atom2.getType() != type2) continue;
                if (test != null && !test.test(atom1, atom2)) continue;
                dr.Ev1Mv2(atom2.getPosition(),atom1.getPosition());
                boundary.nearestImage(dr);
                double r2 = dr.squared();       //compute pair separation
                if(r2 < xMaxSquared) {
                    int index = xDataSource.getIndex(Math.sqrt(r2));  //determine histogram index
                    gSum[index]++;                        //add once for each atom
                }
            }
        }
        callCount++;
    }

    /**
	 * Returns the RDF, averaged over the calls to actionPerformed since the
     * meter was reset or had some parameter changed (xMax or # of bins).
	 */
	public IData getData() {

        if (singleSample) {
            zeroData();
            actionPerformed();
        }
        else {
            if (rData != xDataSource.getData() ||
                    data.getLength() != rData.getLength() ||
                    xDataSource.getXMax() != xMax) {
                reset();
                //that zeroed everything.  just return the zeros.
                return data;
            }
        }
        final double[] y = data.getData();
        long numAtomPairs = 0;
        // we ignore test here because we want the density of atoms pairs that pass
        // the test relative to all atom pairs with nominal density
        if (type1 == null) {
            long numAtoms = box.getLeafList().size();
            numAtomPairs = numAtoms*(numAtoms-1)/2;
        }
        else {
            IAtomList atoms = box.getLeafList();
            long nAtoms1 = 0, nAtoms2 = 0;
            for (int i=0; i<atoms.size(); i++) {
                AtomType t = atoms.get(i).getType();
                if (t == type1) nAtoms1++;
                if (t == type2) nAtoms2++;
            }
            if (type2 == type1) {
                numAtomPairs = nAtoms1 * (nAtoms1 - 1) / 2;
            }
            else {
                numAtomPairs = nAtoms1*nAtoms2;
            }
        }
	    double norm = numAtomPairs * callCount / box.getBoundary().volume();
	    double[] r = rData.getData();
	    double dx2 = 0.5*(xMax - xDataSource.getXMin())/r.length;
	    for(int i=0;i<r.length; i++) {
	        double vShell = space.sphereVolume(r[i]+dx2)-space.sphereVolume(r[i]-dx2);
	        y[i] = gSum[i] / (norm*vShell);
	    }
	    return data;
	}

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray)xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray)xDataSource.getDataInfo();
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    public int getIndependentArrayDimension() {
        return 1;
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
}
