/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.mappedDensity.positionOrientation;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.*;
import etomica.data.types.DataDoubleArray;
import etomica.data.types.DataFunction;
import etomica.molecule.IMolecule;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Measures histogram of molecule positions in z direction, regardless of orientation
 * using Rotenberg's force-averaging formulas (Coles et al (2021))
 * Does not exploit symmetry about centerline between walls
 */
public class MeterDensityProfileForceSimple implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataTag tag;
    protected DataSourceUniform perpSource;
    protected final Box box;
    protected final int d;
    protected double S;
    protected final int nBinsPerp;
    protected final PotentialCompute potentialCompute;
    protected final double T; //temperature

    public MeterDensityProfileForceSimple(PotentialCompute potentialCompute, Box box, double T, int d, int nBinsPerp) {
        this.potentialCompute = potentialCompute;
        this.box = box;
        this.T = T;
        this.d = d;
        this.nBinsPerp = nBinsPerp;
        double L = box.getBoundary().getBoxSize().getX(d);
        perpSource = new DataSourceUniform("perp", Length.DIMENSION, nBinsPerp, 0, L, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        data = new DataFunction(new int[]{nBinsPerp});
        dataInfo = new DataFunction.DataInfoFunction("histogram", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
        S = box.getBoundary().getBoxSize().getX(0);  //area
        if(d == 3) S *= box.getBoundary().getBoxSize().getX(1);
    }

    @Override
    public IData getData() {
        int nMolecules = box.getMoleculeList().size();
        data.E(0);
        double[] y = data.getData();
        double Lz = box.getBoundary().getBoxSize().getX(d);
        potentialCompute.computeAll(true);
        Vector[] forces = potentialCompute.getForces();
        int imol = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            Vector dr = box.getSpace().makeVector();
            IAtomList atoms = molecule.getChildList();
            dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
            box.getBoundary().nearestImage(dr);
            Vector xyz = box.getSpace().makeVector();
            xyz.E(atoms.get(0).getPosition());
            xyz.PEa1Tv1(0.5, dr);
            double zi = xyz.getX(d) + Lz/2;

            imol++;
            Vector f0 = forces[2*imol-2];
            Vector f1 = forces[2*imol-1];
            double fz = f0.getX(d) + f1.getX(d); //total force on molecule in z direction
            fz += wallForce(atoms.get(1).getPosition().getX(d)) + wallForce(atoms.get(0).getPosition().getX(d));

            for(int iz=0; iz<nBinsPerp; iz++) {
                double z = perpSource.getData().getValue(iz);
                if (zi < z) y[iz] += fz / (T * S);
            }
        }
        return data;
    }

    private double wallForce(double zi) {
        double sigma = 1;
        double epsilon = 1;
        double Lz = box.getBoundary().getBoxSize().getX(d);
        double dz = (0.5*Lz - Math.abs(zi));
        if(dz > Math.pow(2, 1. / 6.)) return 0;
        double z3 = 1/(dz*dz*dz);
        double z6 = z3*z3;
        double fw = 48*z6*(z6 - 0.5)/dz;
        return (zi > 0) ? -fw : fw;
    }

    @Override
    public DataTag getTag() {
        return tag;
    }

    @Override
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    @Override
    public DataDoubleArray getIndependentData(int i) {
        return (DataDoubleArray) perpSource.getData();
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataDoubleArray.DataInfoDoubleArray) perpSource.getDataInfo();
    }

    @Override
    public int getIndependentArrayDimension() {
        return 1;
    }

    @Override
    public DataTag getIndependentTag() {
        return null;
    }
}
