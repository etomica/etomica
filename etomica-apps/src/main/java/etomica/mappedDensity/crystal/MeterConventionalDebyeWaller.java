/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.crystal;

import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.normalmode.CoordinateDefinition;
import etomica.space.Vector;
import etomica.units.dimensions.*;

public class MeterConventionalDebyeWaller implements IDataSource, DataSourceIndependent {


    protected final Box box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    protected IDataInfo dataInfo;
     protected Vector rivector;
    protected double msd;
    protected double[] qvector;
    protected int numAtoms;

    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    /**
     * Meter that defines the property being profiled.
     */
    protected final DataTag tag;
    protected CoordinateDefinition latticesite;

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterConventionalDebyeWaller(int numAtoms,double[] qvector, double msd, Box box, CoordinateDefinition latticesite) {
        this.box = box;
        this.qvector = qvector;
        this.numAtoms=numAtoms;
        this.msd = msd;
         this.latticesite = latticesite;
        this.rivector = box.getSpace().makeVector();
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {

        data.E(0);
        double[] y = data.getData();
        IAtomList atoms = box.getLeafList();

        for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
             double qdotrj =0.0;

            for (int i=0;i<3;i++){
                qdotrj =qdotrj +(qvector[i]*rivector.getX(i));         //  double qdotrcap = qvector.dot(rivector)/ri;
            }
     //       System.out.println(qdotrjcap);

            y[0] = y[0]+ qdotrj*qdotrj ;
        }
        y[0]=y[0]/numAtoms;
        return data;
    }

    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) xDataSource.getData();
    }

    public DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataInfoDoubleArray) xDataSource.getDataInfo();
    }

    public int getIndependentArrayDimension() {
        return 1;
    }

    public DataTag getIndependentTag() {
        return xDataSource.getTag();
    }

    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    public void reset() {
        if (box == null) return;
        xDataSource.setXMin(0);
        data = new DataFunction(new int[]{xDataSource.getNValues()});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

}
