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

import java.util.Random;

public class MeterHistogramOrientationMA implements IDataSource, DataSourceIndependent {

    protected DataFunction data;
    protected DataFunction.DataInfoFunction dataInfo;
    protected DataTag tag;
    protected DataSourceUniform perpSource, cosSource;
    protected final Box box;
    protected int d;
    protected int nBinsCos, nBinsPerp;
    protected double S; //area of surface
    protected final double c;
    protected final double[] coshz, cosh2z;
    protected final PotentialCompute potentialCompute;
    protected final double T;

    public MeterHistogramOrientationMA(PotentialCompute potentialCompute, Box box, double T, int d, int nBinsPerp, int nBinsCos) {
        this.potentialCompute = potentialCompute;
        this.box = box;
        this.T = T;
        this.d = d;
        this.nBinsCos = nBinsCos;
        this.nBinsPerp = nBinsPerp;
        double L = box.getBoundary().getBoxSize().getX(d);
        S = box.getBoundary().getBoxSize().getX(0);
        if(d == 3) {
            S *= box.getBoundary().getBoxSize().getX(1);
            //c = ?
            throw new RuntimeException("need to figure this out for 3D");
        } else { //d == 2
            c = 1 / (Math.PI * S * T);
        }
        int D = box.getSpace().D();
        perpSource = new DataSourceUniform("perp", Length.DIMENSION, nBinsPerp, 0, L, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        if(D==2) cosSource = new DataSourceUniform("angle", Null.DIMENSION, nBinsCos, -Math.PI, Math.PI, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        else cosSource = new DataSourceUniform("cosine", Null.DIMENSION, nBinsCos, 0, 1, DataSourceUniform.LimitType.HALF_STEP, DataSourceUniform.LimitType.HALF_STEP);
        coshz = new double[nBinsPerp];
        cosh2z = new double[nBinsPerp];
        for(int i=0; i<nBinsPerp; i++) {
            coshz[i] = Math.cosh(perpSource.getData().getValue(i));
            cosh2z[i] = Math.cosh(2*perpSource.getData().getValue(i));
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
//            System.out.println(fz+", "+torq);

            dr.normalize();
            double drd;//= Math.abs(dr.getX(d)); // z coordinate of orientation
            if (D == 2) drd = Math.atan2(dr.getX(1),dr.getX(0));//angle between -Pi and Pi, with zero on +x axis
            else throw new RuntimeException("Not set up for D != 2");
//            double zi = (0.5*Lz - Math.abs(xyz.getX(d))); // distance from nearer wall
            double zi = xyz.getX(d) + Lz/2;//distance from bottom wall

            // get force and torque on molecule

            //sum over tabulated positions and orientations
            for(int nt = 0; nt<nBinsCos; nt++) {
                double t = cosSource.getData().getValue(nt);
                double ti = drd - t;
                double costi = Math.cos(ti);
                double sinti = Math.sin(ti);
                double cos2ti = costi*costi - sinti*sinti;
                for (int nz = 0; nz < nBinsPerp; nz++) {
                    double z = perpSource.getData().getValue(nz);
                    double coshzi = Math.cosh(zi);
                    double sinhzi = Math.sinh(zi);
                    double cosh2zi = coshzi*coshzi + sinhzi*sinhzi;
                    double sinh2zi = 2*coshzi*sinhzi;
                    double denom = 1 + cos2ti + cosh2zi + cosh2z[nz] - 4*costi*coshz[nz]*coshzi;
                    double tdot = c * 2 * (coshz[nz]*coshzi*sinti - costi*sinti)/denom;
                    double zdot = c * 0.5 * ( -zi/Lz + 2 * (-2*costi*coshz[nz]*sinhzi + sinh2zi)/denom);
                    if(zi > z) zdot -= 0.5 * c; //Heaviside term
                    for(int n = 1; n < 1; n++) {//adjust upper bound to include more terms, or skip entirely
                        tdot += c * Math.cosh(n * z)*(1./Math.tanh(n * Lz) - 1) *
                                    2.*Math.cosh(n * zi) * Math.sin(n * ti);
                    }
                    y[nz*nBinsCos + nt] += -(zdot*fz + tdot*torq);
//                    if(Math.abs(Math.abs(Lz/2-z)-1)<1e-3 && Math.abs(Math.abs(t)-1.0995574287564276)<1e-5) System.out.println(Math.signum(Lz/2-z)+", "+t+", "+ti+", "+tdot);
                }
            }
            //y[perpIdx* cosSource.getNValues() + cosIdx] += inc;
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
}
