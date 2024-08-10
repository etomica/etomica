/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataTag;
import etomica.data.IData;
import etomica.data.IDataInfo;
import etomica.data.IDataSource;
import etomica.data.types.DataDoubleArray;
import etomica.molecule.CenterOfMass;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMasterBonding;
import etomica.potential.compute.PotentialCallback;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMAReal2FD implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1;
    protected double beta_0, dbeta, beta_p, beta_m, omega2;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int nShifts = 0;
    protected int dim;
    protected Vector[] latticePositions;
    protected int numAtoms, nBeads;
    protected Box box;
    protected double EnShift;
    protected Vector[] rOrig;
    protected Vector drShift;
    protected double hbar, mass;
    protected double[] ki_0, gamma_0, dGamma_0, f11_0, f1N_0;
    protected double[][] params_0, params_m, params_p;

    public MeterPIHMAReal2FD(PotentialMasterBonding pmBonding, PotentialCompute pcP1, int nBeads, double omega2, Box box, double temperature, double hbar, double dbeta, int nShifts) {
        int nData = 2;
        this.box = box;
        this.mass = box.getMoleculeList().get(0).getType().getMass(); //RP mass
        this.nShifts = nShifts;
        this.hbar = hbar;
        this.nBeads = nBeads;
        this.dbeta = dbeta;
        this.beta_0 = 1/temperature;
        this.beta_p = beta_0 + dbeta;
        this.beta_m = beta_0 - dbeta;
        this.omega2 = omega2;

        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI", Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        numAtoms = box.getMoleculeList().size();
        dim = box.getSpace().D();
        rOrig = box.getSpace().makeVectorArray(box.getLeafList().size());
        latticePositions = box.getSpace().makeVectorArray(numAtoms);
        for (int i = 0; i < latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

        params_0 = computeMappingParams(beta_0);
        params_m = computeMappingParams(beta_m);
        params_p = computeMappingParams(beta_p);

        ki_0     = params_0[0];
        gamma_0  = params_0[1];
        dGamma_0 = params_0[2];
        f11_0    = params_0[3];
        f1N_0    = params_0[4];

        double En_ho_stage = dim*nBeads/2.0/beta_0;
        for (int k = 0; k < nBeads; k++) {
            En_ho_stage -= dim*gamma_0[k];
        }
        //COM
        En_ho_stage -= dim/2.0/beta_0;
        System.out.println(" En_ho_stage:  " + En_ho_stage);

        this.EnShift = 0;
        this.drShift = box.getSpace().makeVector();

        if (nShifts < 0) throw new RuntimeException("nShifts cannot be negative");
        int ns = nShifts + 1;
        if (ns > nBeads || (nBeads/ns)*ns != nBeads) {
            throw new RuntimeException("# of beads must be a multiple of (# of shifts + 1)");
        }
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        for (int i = 0; i < x.length; i++) x[i] = 0;

        if (box.getMoleculeList().size() > 1) drShift = computeShift();

        for (IMolecule molecule : box.getMoleculeList()) {
            for (IAtom atom : molecule.getChildList()) {
                rOrig[atom.getLeafIndex()].E(atom.getPosition());
            }
        }
        double En = dim*numAtoms*nBeads/2.0/beta_0;
        double Cvn = dim*numAtoms*nBeads/2.0/beta_0/beta_0;
        for (int i = 0; i < nBeads; i++) {
            En -= dim*numAtoms*gamma_0[i];
            Cvn += dim*numAtoms*dGamma_0[i];
        }
        if (numAtoms > 1 && omega2 != 0) { //COM
            En -= dim/2.0/beta_0;
            Cvn -= dim/2.0/beta_0/beta_0;
        }

        double dbUdb   = getdbUdb(beta_0, params_0);
        double dbUdb_p = getdbUdb(beta_p, params_p);
        double dbUdb_m = getdbUdb(beta_m, params_m);
        double d2bUdb2 = (dbUdb_p - dbUdb_m)/(2*dbeta);

        x[0] = En + dbUdb - EnShift;
        x[1] = Cvn - d2bUdb2  + x[0]*x[0];

        return data;
    }

    private double getdbUdb(double beta, double[][] params) {
        if (beta != beta_0) scaleCoord(beta, params);
        double[] gamma = params[1];
        double[] f11   = params[3];
        double[] f1N   = params[4];
        double[] df11  = params[5];
        double[] df1N  = params[6];

        pcP1.computeAll(true);
        double U = pcP1.getLastEnergy();
        Vector[] forcesU = pcP1.getForces();

        pmBonding.computeAll(true);
        double facK = Math.pow(beta_0/beta, 2);
        double K = facK*pmBonding.getLastEnergy();

        Vector[] forcesK = box.getSpace().makeVectorArray(box.getLeafList().size());
        for (int i = 0; i < forcesK.length; i++){
            forcesK[i].Ea1Tv1(facK, pmBonding.getForces()[i]);
        }

        Vector tmp_r = box.getSpace().makeVector();
        Vector drj = box.getSpace().makeVector();
        Vector v = box.getSpace().makeVector();
        Vector v0 = box.getSpace().makeVector();
        Vector vPrev = box.getSpace().makeVector();
        Vector dr0 = box.getSpace().makeVector();

        IMoleculeList molecules = box.getMoleculeList();

        double dbUdb = 0;
        int ns = nShifts + 1;
        for (int indexShift = 0; indexShift < nBeads; indexShift += nBeads/ns) {
            for (int i = 0; i < molecules.size(); i++) {
                IAtomList beads = molecules.get(i).getChildList();
                dr0.Ev1Mv2(beads.get(indexShift).getPosition(), latticePositions[i]);
                dr0.PE(drShift);
                Vector drPrev = box.getSpace().makeVector();
                for (int j = 0; j < beads.size(); j++) {
                    int aj = (j + indexShift) % beads.size();
                    IAtom atomj = beads.get(aj);
                    drj.Ev1Mv2(atomj.getPosition(), latticePositions[i]);
                    drj.PE(drShift);
                    tmp_r.E(drj);
                    if (j > 0) {
                        tmp_r.PEa1Tv1(-f11[j], drPrev);
                        tmp_r.PEa1Tv1(-f1N[j], dr0);
                    }
                    v.Ea1Tv1(gamma[j], tmp_r);
                    if (j > 0) {
                        v.PEa1Tv1(df11[j], drPrev);
                        v.PEa1Tv1(df1N[j], dr0);
                        v.PEa1Tv1(f11[j], vPrev);
                        v.PEa1Tv1(f1N[j], v0);
                    }
                    int jj = atomj.getLeafIndex();
                    dbUdb -= beta*forcesU[jj].dot(v) + beta*forcesK[jj].dot(v);
                    drPrev.E(drj);
                    if (j == 0) {
                        v0.E(v);
                    }
                    vPrev.E(v);
                } // beads
            }//mol
        }//shifts
        dbUdb /= ns;
        dbUdb += U - K;

        // If atoms scaled, bring them back to their original coords
        if (beta != beta_0) {
            for (IMolecule molecule : box.getMoleculeList()) {
                for (IAtom atom : molecule.getChildList()) {
                    atom.getPosition().E(rOrig[atom.getLeafIndex()]);
                }
            }
        }

        return dbUdb;
    }

    private void scaleCoord(double beta, double[][] params) {
        double[] ki  = params[0];
        double[] f11 = params[3];
        double[] f1N = params[4];

        Vector dr0 = box.getSpace().makeVector();
        Vector dri = box.getSpace().makeVector();
        Vector drPrev0 = box.getSpace().makeVector();
        Vector drPrev = box.getSpace().makeVector();
        Vector[] u = box.getSpace().makeVectorArray(nBeads);

        for (IMolecule molecule : box.getMoleculeList()) {
            IAtomList atoms = molecule.getChildList();
            dr0.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[molecule.getIndex()]);
            dr0.PE(drShift);
            for (IAtom atom : atoms) {
                // r->u
                int i = atom.getIndex();
                dri.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                dri.PE(drShift);
                u[i].E(dri);
                if (i > 0) {
                    u[i].PEa1Tv1(-f1N_0[i], dr0);
                    u[i].PEa1Tv1(-f11_0[i], drPrev0);
                }
                drPrev0.E(dri);

                //scaling: u->u'
                u[i].TE(Math.sqrt(beta_0*ki_0[i]/(beta*ki[i])));

                //u'->r'
                atom.getPosition().E(latticePositions[molecule.getIndex()]);
                atom.getPosition().PE(u[i]);

                if (i > 0) {
                    atom.getPosition().PEa1Tv1(f1N[i], u[0]);
                    atom.getPosition().PEa1Tv1(f11[i], drPrev);
                }
                drPrev.Ev1Mv2(atom.getPosition(), latticePositions[molecule.getIndex()]);
                atom.getPosition().ME(drShift);
            } // a
        } // m
    }

    private double[][] computeMappingParams(double beta) {
        double[] ki = new double[nBeads];
        double[] gamma = new double[nBeads];
        double[] dGamma = new double[nBeads];
        double[] f11 = new double[nBeads];
        double[] f1N = new double[nBeads];
        double[] df11 = new double[nBeads];
        double[] df1N = new double[nBeads];

        double omegaN = nBeads/(hbar*beta);
        double omegaN2 = omegaN*omegaN;
        double D = 2 + omega2/omegaN2;
        double alpha = Math.log(D/2 + Math.sqrt(D*D/4 - 1));
        double dAlpha = 2.0/beta/Math.sqrt(1.0 + 4.0*omegaN2/omega2);
        double dAlpha2 = dAlpha*dAlpha;
        double d2Alpha = -1.0/4.0*beta*dAlpha*dAlpha*dAlpha;
        double sinhA = Math.sinh(alpha);
        double coshA = Math.cosh(alpha);
        double sinhNA = Math.sinh(nBeads*alpha);
        double coshhNA = Math.cosh(nBeads*alpha);

        ki[0] = 2*mass/nBeads*omegaN2*sinhA*Math.tanh(nBeads*alpha/2.0);
        gamma[0] = 1.0/2.0/beta - dAlpha/2.0*(coshA/sinhA + nBeads/sinhNA);
        dGamma[0] = -1.0/2.0/beta/beta - 1.0/2.0*d2Alpha*(coshA/sinhA + nBeads/sinhNA)
                  + 0.5*dAlpha2*(1/sinhA/sinhA + nBeads*nBeads/sinhNA*coshhNA/sinhNA);

        for (int i = 1; i < nBeads; i++) {
            double sinhNmiA = Math.sinh((nBeads-i)*alpha);
            double coshNmiA = Math.cosh((nBeads-i)*alpha);
            double sinhNmip1A = Math.sinh((nBeads-i+1)*alpha);
            double coshNmip1A = Math.cosh((nBeads-i+1)*alpha);

            double sinhRatio = sinhNmip1A/sinhNmiA;
            ki[i] = mass/nBeads*omegaN2*sinhRatio;
            gamma[i] = 1.0/2.0/beta - dAlpha/2.0*(coshNmip1A/sinhNmip1A - (nBeads-i)*sinhA/sinhNmip1A/sinhNmiA);
            dGamma[i] = -1.0/2.0/beta/beta-d2Alpha/2.0*coshNmip1A/sinhNmip1A + (nBeads-i)/2.0*d2Alpha*sinhA/sinhNmip1A/sinhNmiA
                    + dAlpha2/2.0*(nBeads-i+1)/sinhNmip1A/sinhNmip1A
                    + dAlpha2/2.0*(nBeads-i)*coshA/sinhNmip1A/sinhNmiA
                   -dAlpha2/2.0*(nBeads-i)*(nBeads-i+1)*sinhA/sinhNmip1A*coshNmip1A/sinhNmip1A/sinhNmiA
                   -dAlpha2/2.0*(nBeads-i)*(nBeads-i)*sinhA/sinhNmip1A/sinhNmiA*coshNmiA/sinhNmiA;

            f11[i] = sinhNmiA/sinhNmip1A;
            f1N[i] = sinhA/sinhNmip1A;
            df11[i] = dAlpha/sinhNmip1A*((nBeads-i)*coshNmiA-(nBeads-i+1)*sinhNmiA*coshNmip1A/sinhNmip1A);
            df1N[i] = dAlpha/sinhNmip1A*(coshA-(nBeads-i+1)*sinhA*coshNmip1A/sinhNmip1A);
        }

        return new double[][] {ki, gamma, dGamma, f11, f1N, df11, df1N};
    }

    public void setEnShift(double E) {
        this.EnShift = E;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    protected Vector computeShift() {
        if (box.getMoleculeList().size() == 1) return box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        Vector drShift = box.getSpace().makeVector();
        drShift.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[0]);
        drShift.TE(-1);
        return drShift;
    }
}