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
import etomica.potential.IteratorDirective;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialMaster;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.*;

/**
 *
 * Meter for calculating force sampling singlet density of fluids
 *
 */
public class MeterProfileForceSum implements IDataSource, DataSourceIndependent, AtomLeafAgentManager.AgentSource<Vector> {

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

    /**
     * Default constructor sets profile along the y-axis, with 100 histogram points.
     */
    public MeterProfileForceSum(Box box, PotentialMaster potentialMaster, double temperature) {
        this.box = box;
        this.temperature = temperature;
        Space space = box.getSpace();
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

    /**
     * Returns the profile for the current configuration.
     */
    public IData getData() {
        pc.reset();
        potentialMaster.calculate(box, id, pc);
        Boundary boundary = box.getBoundary();
        data.E(0);
        double[] y = data.getData();
        IAtomList atoms = box.getLeafList();
        double L = box.getBoundary().getBoxSize().getX(profileDim);
        double dz = L / xDataSource.getNValues();
        for (IAtom atom : atoms) {
            double fz = agentManager.getAgent(atom).getX(profileDim);
            double zi = atom.getPosition().getX(profileDim);
            if(zi<-L/2+dz/2){continue;}
            double binzi = zi + 0.5 * dz;
            if (binzi > L / 2) continue;
            int izi = xDataSource.getIndex(binzi);
            for (int i = izi; i < y.length; i++) {
                y[i] += fz;
            }
        }

        double area = box.getBoundary().volume() / L;
        data.TE(1 / (temperature * area));
        double sumRho = 0;
        for (int i = 0; i < y.length; i++) {
            sumRho += y[i];
        }
        double rhoCalc = sumRho / y.length;
        double drho = atoms.size() / box.getBoundary().volume() - rhoCalc;
        data.PE(drho);
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
        dataInfo = new DataInfoFunction("Force Sum Profile", new CompoundDimension(new Dimension[]{Quantity.DIMENSION, Volume.DIMENSION}, new double[]{1, -1}), this);
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
