/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.mappedDensity;

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
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.units.dimensions.*;

public class MeterProfileMappedAvg implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

    protected final PotentialMaster potentialMaster;
    protected final IteratorDirective id;
    protected final PotentialCalculationForceSum pc;
    protected final AtomLeafAgentManager<Vector> agentManager;
    protected final Box box;
    protected DataSourceUniform xDataSource;
    protected DataFunction data;
    protected IDataInfo dataInfo;
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

    protected final FunctionDifferentiable c;
    protected Behavior behavior;
    protected double zidotz;

    public enum Behavior {
        NORMAL, P, ZIDOT, DZIDOT
    }

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfileMappedAvg(Box box, PotentialMaster potentialMaster, double temperature, FunctionDifferentiable c) {
        this.box = box;
        this.temperature = temperature;
        this.c = c;
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

    public double zidot(double z, double zi) {
        double L = box.getBoundary().getBoxSize().getX(profileDim);
        double pzi = c.df(1, zi);
        double pz = c.df(1, z);
        double q = c.f(L / 2);
        double cz = c.f(z);
        double czi = c.f(zi);
  //      double czplusLby2 = c.f(z+L/2);
  //      double czminusLby2 = c.f(z-L/2);

        // c(zi) computed starting from z
        // p(zi) zidot(zi) = (p(z) zidot(z+)) - beta p(z) c(zi)/q
        // zi=z- => p(z) zidot(z-) = (p(z) zidot(z+)) - beta p(z)
        // z+ - z- => p(z) zidot(z+) - p(z) zidot(z-) = + beta p(z)
        // zidot(z+) - zidot(z-) = beta
        // zidot(z+) = beta/2
        //
        // p(zi) zidot(zi) = beta p(z) / 2 - beta p(z) c(zi)/q
        //                 = beta p(z) (1/2 - c(zi)/q)
        // zidot(zi) = beta p(z)/p(zi) (1/2 - c(zi)/q)

        // our c(zi) wasn't actually computed starting from z
//////////////////////////
        double x = (zi > z) ? (czi - cz) / q : ((czi - cz) / q + 1);
        return pz / pzi * (0.5 - x) / temperature;

 //       double pminusLby2= c.df(1, -L/2);
 //       double zidotminusLby2=-0.0725;
   //     zidotminusLby2=pz*(cz/q-1/2)/(pminusLby2*temperature);
  //      if (zi >= z) {zidotminusLby2=-pz*(1-(czplusLby2/q))/(pminusLby2*temperature);} else {zidotminusLby2=pz*((czminusLby2/q))/(pminusLby2*temperature);}
 //       double heavisidei;
 //       if (zi >= z) {heavisidei=1;} else {heavisidei=0;}
  //      return ((pminusLby2*zidotminusLby2/pzi)+(pz*((heavisidei)-(czi/q))/(temperature*pzi)));
   //     return ((((heavisidei)-(zi/L))/(temperature)));

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
        double dz = L / xDataSource.getNValues();
        if (behavior == Behavior.ZIDOT) {
            for (int i = 0; i < y.length; i++) {
                double zi = -L / 2 + i * dz;
                if (Math.abs(zidotz - zi) < dz * 0.01) {
                    y[i] = Double.NaN;
                    continue;
                }
                y[i] = zidot(zidotz, zi);
            }
            return data;
        }
        if (behavior == Behavior.DZIDOT) {
            for (int i = 0; i < y.length; i++) {
                double zi = -L / 2 + i * dz;
                if (Math.abs(zidotz - zi) < dz * 0.01) {
                    y[i] = Double.NaN;
                    continue;
                }
                y[i] = -(c.df(1, zi + dz * 0.01) * zidot(zidotz, zi + dz * 0.01) - c.df(1, zi - dz * 0.01) * zidot(zidotz, zi - dz * 0.01)) / (0.02 * dz) / c.df(1, zi);
            }
            return data;
        }
        if (behavior != Behavior.P) {
            for (IAtom atom : atoms) {
                double fz = agentManager.getAgent(atom).getX(profileDim);
                double zi = atom.getPosition().getX(profileDim);
                for (int i = 0; i < y.length; i++) {
                    double z = -L / 2 + (i + 0.5) * dz;
                    // dphi/dz = -T d(ln(p))/dz
                    //         = -T / p dp/dz
                    // - T p dp/dz
                  y[i] -= (fz - c.df(2, zi) / c.df(1, zi) *temperature) * zidot(z, zi);

                }
            }
        }
//        data.TE(0);

        int N = atoms.size();
        double q = c.f(L / 2);
        for (int i = 0; i < y.length; i++) {
            double z = -L / 2 + (i + 0.5) * dz;
            y[i] += N * c.df(1, z) / q;
        }
        double area = box.getBoundary().volume() / L;
        data.TE(1 / area);
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
        double halfBox = 0.5 * boundary.getBoxSize().getX(profileDim);
        xDataSource.setXMin(-halfBox);
        xDataSource.setXMax(halfBox);

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
