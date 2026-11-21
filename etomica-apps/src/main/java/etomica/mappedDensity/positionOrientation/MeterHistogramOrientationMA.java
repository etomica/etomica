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
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

public class MeterHistogramOrientationMA implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataTag tag;
    protected IDataSource perpSource;
    protected DataSourceUniform cosSource;
    protected final Box box;
    protected int d;
    protected int nBinsCos, nBinsPerp;
    protected double S; //area of surface
    protected final double c;
    protected final double[][] coshz, cosh2z, coshzM, cosh2zM;
    protected final PotentialCompute potentialCompute;
    protected final double T;
    protected final double[] perps;
    protected final int nmax = 2;
    protected final double[] corrn = new double[nmax];

    public MeterHistogramOrientationMA(PotentialCompute potentialCompute, Box box, double T, int d, int nBinsCos, double[] perps) {
        this(potentialCompute, box, T, d, nBinsCos, perps, perps.length);
    }

    public MeterHistogramOrientationMA(PotentialCompute potentialCompute, Box box, double T, int d, int nBinsPerp, int nBinsCos) {
        this(potentialCompute, box, T, d, nBinsCos, null, nBinsPerp);
    }

    protected MeterHistogramOrientationMA(PotentialCompute potentialCompute, Box box, double T, int d, int nBinsCos, double[] perps, int nBinsPerp) {
        this.potentialCompute = potentialCompute;
        this.box = box;
        this.T = T;
        this.d = d;
        this.nBinsCos = nBinsCos;
        this.nBinsPerp = nBinsPerp;
        this.perps = perps;
        double L = box.getBoundary().getBoxSize().getX(d);
        S = box.getBoundary().getBoxSize().getX(0);
        if(d == 3) {
            S *= box.getBoundary().getBoxSize().getX(1);
            //c = ?
            throw new RuntimeException("need to figure this out for 3D");
        } else { //d == 2
            c = 1 / (Math.PI * S * L * T);
        }
        int D = box.getSpace().D();
        if (perps == null) {
            perpSource = new DataSourceUniform("perp", Length.DIMENSION, nBinsPerp, 0, L/2, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        }
        else {
            perpSource = new DataSourceStatic("perp", Length.DIMENSION, perps);
        }
        if(D==2) cosSource = new DataSourceUniform("angle", Null.DIMENSION, nBinsCos, -Math.PI, Math.PI, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        else cosSource = new DataSourceUniform("cosine", Null.DIMENSION, nBinsCos, 0, 1, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        coshz = new double[nBinsPerp][2];
        cosh2z = new double[nBinsPerp][2];
        coshzM = new double[nBinsPerp][2];
        cosh2zM = new double[nBinsPerp][2];

        for(int i=0; i<nBinsPerp; i++) {
            for (int j=0; j<2; j++) {
                double z = perpSource.getData().getValue(i);
                if (j == 0) z = L/2 + z;  // from L/2 to L
                else z = L/2 - z;         // from L/2 to 0
                double zred = z/L;
                coshz[i][j] = Math.cosh(zred);
                cosh2z[i][j] = Math.cosh(2 * zred);
                coshzM[i][j] = Math.cosh(2-zred);
                cosh2zM[i][j] = Math.cosh(2 * (2-zred));

            }
        }

        for(int n=0; n<nmax; n++) {
            corrn[n] = 1.0/Math.sinh(n+1) - 2*Math.exp(-(n+1));
        }

        data = new DataFunction(new int[]{nBinsPerp,nBinsCos});
        dataInfo = new DataFunction.DataInfoFunction("histogram", Null.DIMENSION, this);
        tag = new DataTag();
        dataInfo.addTag(tag);
    }

    @Override
    public IData getData() {
        int nMolecules = box.getMoleculeList().size();

        data.E(nMolecules/(2*Math.PI * box.getBoundary().volume()));
        double[] y = data.getData();
        double Lz = box.getBoundary().getBoxSize().getX(d);
        int D = box.getSpace().D();

        potentialCompute.computeAll(true);
        Vector[] forces = potentialCompute.getForces();
        Vector torque = box.getSpace().makeVector();
        int imol = 0;
        for (IMolecule molecule : box.getMoleculeList()) {
            imol++;

            //get position and direction of molecule
            Vector dr = box.getSpace().makeVector(); //dr will point from atom 0 to atom 1
            IAtomList atoms = molecule.getChildList();
            dr.Ev1Mv2(atoms.get(1).getPosition(), atoms.get(0).getPosition());
            box.getBoundary().nearestImage(dr);
            Vector xyz = box.getSpace().makeVector();
            xyz.E(atoms.get(0).getPosition());
            xyz.PEa1Tv1(0.5, dr); //xyz is location of center of pair

            Vector f0 = forces[2*imol-2];
            Vector f1 = forces[2*imol-1];
            f0.setX(d,f0.getX(d) + wallForce(atoms.get(0).getPosition().getX(d)));
            f1.setX(d,f1.getX(d) + wallForce(atoms.get(1).getPosition().getX(d)));

            double fz = f0.getX(d) + f1.getX(d); //total force on molecule in z direction
            if(fz == 0.0) continue;

            torque.Ev1Mv2(f1,f0); // torque = dr/2 x f1 + (-dr/2) x f0 = dr/2 x (f1 - f0)

            //special for 2D
            if(D == 3) throw new RuntimeException("need to update torque for 3D");
            double torq = 0.5*(torque.getX(1)*dr.getX(0) - torque.getX(0)*dr.getX(1)); //fy*rx - fx*ry

            dr.normalize();
            double drd;//= Math.abs(dr.getX(d)); // z coordinate of orientation
            if (D == 2) drd = Math.atan2(dr.getX(1),dr.getX(0));//angle between -Pi and Pi, with zero on +x axis
            else throw new RuntimeException("Not set up for D != 2");
//            double zi = (0.5*Lz - Math.abs(xyz.getX(d))); // distance from nearer wall
            double zi = xyz.getX(d) + Lz/2;//distance from bottom wall
            double ziRed = zi/Lz;
            double coshzi = Math.cosh(ziRed);
            double sinhzi = Math.sinh(ziRed);
            double cosh2zi = coshzi * coshzi + sinhzi * sinhzi;

            for(int nt = 0; nt<nBinsCos; nt++) {
                double t = cosSource.getData().getValue(nt); //tabulated orientation, angle wrt x axis
                double ti = drd - t; //orientation of molecule relative to tabulated orientation
                double costi = Math.cos(ti);
                double sinti = Math.sin(ti);
                double cos2ti = costi*costi - sinti*sinti;
                for (int nz = 0; nz < nBinsPerp; nz++) {
                    for (int j=0; j<2; j++) {
                        double z = 0;
                        if (nmax > 0) {
                            z = perpSource.getData().getValue(nz);
                            if (j == 0) z = Lz / 2 + z;  // from Lz/2 to Lz
                            else z = Lz / 2 - z;          // from Lz/2 to 0
                        }
                        double denom = 1 + cos2ti + cosh2zi + cosh2z[nz][j] - 4 * costi * coshz[nz][j] * coshzi;
                        double denomM = 1 + cos2ti + cosh2zi + cosh2zM[nz][j] - 4 * costi * coshzM[nz][j] * coshzi;
                        double tdot = 2*sinti * (coshz[nz][j] * coshzi - costi) / denom;//diverges for ri->r, ti->0
                        if (Math.abs(tdot) / Math.PI > Double.MAX_VALUE) {
                            tdot = 0;//truncate if diverging
                            System.out.println("truncating tdot");
                        }
                        tdot += 2*sinti * (coshzM[nz][j] * coshzi - costi) / denomM;//doesn't diverge
                        double zdot = 2*sinhzi * (-costi * coshz[nz][j] + coshzi) / denom;//diverges for ri->r, ti->0
                        if (Math.abs(zdot) / Math.PI > Double.MAX_VALUE) {
                            zdot = 0;//truncate if diverging
                            System.out.println("truncating zdot");
                        }
                        zdot += -ziRed + 2*sinhzi * (-costi * coshzM[nz][j] + coshzi) / denomM;//doesn't diverge
                        for (int n = 1; n<=nmax; n++) {//correction from approximation of 1/sinh
                            double coshn1r = Math.cosh(n*(1-z/Lz));
                            zdot += (-2*Math.sinh(n*ziRed)*coshn1r)*corrn[n-1]*Math.cos(n*ti);
                            tdot += (+2*Math.cosh(n*ziRed)*coshn1r)*corrn[n-1]*Math.sin(n*ti);
                        }
                        // 0.5 multiplier is introduced for extension of [0,Pi] solution to [-Pi,Pi]
                        // another 0.5 multiplier accounts for double counting j = 0, 1, sources
                        y[nz * nBinsCos + nt] += -c * 0.5 * 0.5 * (Lz * zdot * fz + tdot * torq);
                    }
                }
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
        return (DataDoubleArray) (i==0 ? perpSource.getData() : cosSource.getData());
    }

    @Override
    public DataDoubleArray.DataInfoDoubleArray getIndependentDataInfo(int i) {
        return (DataDoubleArray.DataInfoDoubleArray) (i==0 ? perpSource.getDataInfo() : cosSource.getDataInfo());
    }

    @Override
    public int getIndependentArrayDimension() {
        return 2;
    }

    @Override
    public DataTag getIndependentTag() {
        return null;
    }

    public class DataSourceStatic implements IDataSource {
        protected final DataDoubleArray data;
        protected final DataDoubleArray.DataInfoDoubleArray dataInfo;
        protected final DataTag tag;

        public DataSourceStatic(String label, Dimension dim, double[] x) {
            data = new DataDoubleArray(new int[]{x.length}, x);
            tag = new DataTag();
            dataInfo = new DataDoubleArray.DataInfoDoubleArray(label, dim, new int[]{x.length});
            dataInfo.addTag(tag);
        }

        @Override
        public IData getData() {
            return data;
        }

        @Override
        public DataTag getTag() {
            return tag;
        }

        @Override
        public IDataInfo getDataInfo() {
            return dataInfo;
        }
    }
}
