/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.mappedDensityfromlatticesite;

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
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.*;

public class MeterConventionalqdotr implements IDataSource, DataSourceIndependent {


    protected final Box box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected double Rmax;
    protected Vector rivector;
    protected double msd;

    /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    /**
     * Meter that defines the property being profiled.
     */
    protected final DataTag tag;
    protected CoordinateDefinition latticesite;
    protected double[] qvector;

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterConventionalqdotr(double[] qvector,double msd, Box box, CoordinateDefinition latticesite) {
        this.box = box;
        this.qvector = qvector;

        this.msd = msd;
        this.Rmax = Math.sqrt(msd)*6;
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
//sim.coordinateDefinition - latticesites
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     */
    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    //       double pminusLby2= c.df(1, -L/2);
    //       double zidotminusLby2=-0.0725;
    //     zidotminusLby2=pz*(cz/q-1/2)/(pminusLby2*temperature);
    //      if (zi >= z) {zidotminusLby2=-pz*(1-(czplusLby2/q))/(pminusLby2*temperature);} else {zidotminusLby2=pz*((czminusLby2/q))/(pminusLby2*temperature);}
    //      double heavisidei;
    //      if (zi >= z) {heavisidei=1;} else {heavisidei=0;}
    //      return ((pminusLby2*zidotminusLby2/pzi)+(pz*((heavisidei)-(czi/q))/(temperature*pzi)));
    //    return ((((heavisidei)-(zi/L))/(temperature)));


    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {

        data.E(0);
        double[] y = data.getData();
        IAtomList atoms = box.getLeafList();
        double dz = Rmax / xDataSource.getNValues();

        for (IAtom atom : atoms) {
            rivector.Ev1Mv2(atom.getPosition(), latticesite.getLatticePosition(atom));
            box.getBoundary().nearestImage(rivector);
            double ri = Math.sqrt(rivector.squared());
            int i = xDataSource.getIndex(ri);
            y[i] ++;
        }
     for (int bin =0;bin<y.length;bin++ )   {
            double Rbinstart=bin*dz;
         double Rbinend=(bin+1)*dz;
        y[bin] = y[bin]/(4*Math.PI*(Rbinend*Rbinend*Rbinend-Rbinstart*Rbinstart*Rbinstart)/3);
      //   System.out.println("qqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqqq "+q);

     }
//        data.TE(0);
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

        Boundary boundary = box.getBoundary();
        xDataSource.setXMin(0);
        xDataSource.setXMax(Rmax);

        data = new DataFunction(new int[]{xDataSource.getNValues()});
        dataInfo = new DataInfoFunction("Mapped Average Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
        dataInfo.addTag(tag);
    }

    public DataSourceUniform getXDataSource() {
        return xDataSource;
    }

}
