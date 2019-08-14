/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity.crystal;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.DataSourceUniform.LimitType;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataDoubleArray.DataInfoDoubleArray;
import etomica.data.types.DataFunction;
import etomica.data.types.DataFunction.DataInfoFunction;
import etomica.math.function.FunctionDifferentiable;
import etomica.normalmode.CoordinateDefinition;
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.*;

//Meter for calculating 1D mapped averaging averaging sdf with 1D mapping
//probability=gaussian

public class MeterMappedAvg implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected final PotentialCalculationForceSum pc;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final Box box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    protected IDataInfo dataInfo;
    protected double Rmax;
    protected Vector rivector;
     /**
     * Vector describing the orientation of the profile.
     * For example, (1,0) is along the x-axis.
     */
    protected int profileDim;
    /**
     * Meter that defines the property being profiled.
     */
    protected final DataTag tag;
    protected double temperature;
protected CoordinateDefinition latticesite;
    protected final FunctionDifferentiable c;
    protected Behavior behavior;
    protected double zidotz;
    protected double msd;

    public enum Behavior {
        NORMAL, P, ZIDOT, DZIDOT
    }

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterMappedAvg(double msd,Box box, PotentialMaster potentialMaster, double temperature, FunctionDifferentiable c, CoordinateDefinition latticesite) {
        this.box = box;
        this.temperature = temperature;
        this.c = c;
        this.msd = msd;
        this.Rmax = Math.sqrt(msd)*6;
this.latticesite=latticesite;
this.rivector =box.getSpace().makeVector();
        xDataSource = new DataSourceUniform("x", Length.DIMENSION);
        tag = new DataTag();
        xDataSource.setTypeMax(LimitType.HALF_STEP);
        xDataSource.setTypeMin(LimitType.HALF_STEP);
        this.potentialMaster = potentialMaster;
        id = new IteratorDirective();
        pc = new PotentialCalculationForceSum();
        agentManager = new AtomLeafAgentManager<Vector>(this, box);
        pc.setAgentManager(agentManager);
    }

    public void setBehavior(Behavior b) {
        behavior = b;
    }

    public void setZidotZ(double z) {
        zidotz = z;
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
    public int getProfileDim() {
        return profileDim;
    }

    /**
     * Accessor method for vector describing the direction along which the profile is measured.
     * Each atom position is dotted along this vector to obtain its profile abscissa value.
     * The given vector is converted to a unit vector, if not already.
     */
    public void setProfileDim(int dim) {
        profileDim = dim;
        reset();
    }

    public double ridot(double r, double ri) {
        double pri = c.df(1, ri);
        double pr = c.df(1, r);
        double q = 1;
        double cri = c.f(ri);
//////////////////////////
        double heavisidei;
        if (ri >= r) {heavisidei=1;} else {heavisidei=0;}

    //    System.out.println(" r: " + r+" ri: " + ri+" rsqp*ridot: " +  (r*r*pr*pr* ((heavisidei/(4*Math.PI)) - cri/q)) / (temperature*ri*ri*pri) );

        return (pr  * ((heavisidei/(4*Math.PI)) - cri/q)) / (temperature*ri*ri*pri);

    }

    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        data.E(0);
        double[] y = data.getData();
        IAtomList atoms = box.getLeafList();
        double L = box.getBoundary().getBoxSize().getX(profileDim);
        double dz = Rmax / xDataSource.getNValues();
        if (behavior == Behavior.ZIDOT) {
            for (int i = 0; i < y.length; i++) {
                double zi = i * dz;
                if (Math.abs(zidotz - zi) < dz * 0.01) {
                    y[i] = Double.NaN;
                    continue;
                }

                y[i] = ridot(zidotz, zi);
                System.out.println(" i "+i+" y[i] "+y[i]);

            }
            return data;
        }
        if (behavior == Behavior.DZIDOT) {
            for (int i = 0; i < y.length; i++) {
                double zi =  i * dz;
                if (Math.abs(zidotz - zi) < dz * 0.01) {
                    y[i] = Double.NaN;
                    continue;
                }
                y[i] = -(c.df(1, zi + dz * 0.01) * ridot(zidotz, zi + dz * 0.01) - c.df(1, zi - dz * 0.01) * ridot(zidotz, zi - dz * 0.01)) / (0.02 * dz) / c.df(1, zi);

                System.out.println(" i "+i+" y[i] "+y[i]);

            }

            return data;
        }
        if (behavior != Behavior.P) {
            for (IAtom atom : atoms) {
                 rivector.Ev1Mv2(atom.getPosition(),latticesite.getLatticePosition(atom)) ;
                 box.getBoundary().nearestImage(rivector);
                double ri = Math.sqrt(rivector.squared()) ;
                double fr = agentManager.getAgent(atom).dot(rivector)/ri;

     //           System.out.println(y[0]);


                for (int i = 0; i < y.length; i++) {
                    double r =  (i + 0.5) * dz;

         //      System.out.println(" r "+r+ " ri "+ri+ " ensemble "+(fr - c.df(2, ri) / c.df(1, ri) * temperature) * ridot(r, ri)/(r*r)+ " y[i] "+y[i]);


               //     y[i] -= (fr - (temperature*c.df(2, ri) / c.df(1, ri)) ) * ridot(r, ri);


                }
            }

        }
//        data.TE(0);

        int N = atoms.size();
        double q = 1;
        for (int i = 0; i < y.length; i++) {
            double r =  (i + 0.5) * dz;
            y[i] += N * c.df(1, r) / q;

        }

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

    @Override
    public Vector makeAgent(IAtom a, Box agentBox) {
        return box.getSpace().makeVector();
    }

    @Override
    public void releaseAgent(Vector agent, IAtom atom, Box agentBox) {

    }
}
