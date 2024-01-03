/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.normalmode;

import etomica.atom.IAtom;
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
import etomica.space.Boundary;
import etomica.space.Tensor;
import etomica.space.Vector;
import etomica.units.dimensions.Null;

public class MeterPIHMAReal2Inf implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1harm, pcP1ah;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected int nShifts = 0;
    protected int dim;
    protected Vector[] rdot;
    protected double drdotHdrdot;
    protected Vector[] latticePositions;
    protected int numAtoms, nBeads;
    protected double[] gamma, dGamma, f11, f1N, df11, df1N, d2f11, d2f1N;
    protected Box box;
    protected double hbar, omega, omegaSample, temperature, beta;
    protected double R, sinhR, cothR, massRing;

    public MeterPIHMAReal2Inf(PotentialMasterBonding pmBonding, PotentialCompute pcP1harm,PotentialCompute pcP1ah, double temperature, int nBeads, double omega,double omegaSample, Box box, double hbar) {
        int nData = 2;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);

        this.temperature = temperature;
        this.beta = 1/temperature;
        this.nBeads = nBeads;
        this.hbar = hbar;
        this.omega = omega;
        this.omegaSample = omegaSample;
        this.pmBonding = pmBonding;
        this.pcP1harm = pcP1harm;
        this.pcP1ah = pcP1ah;
        this.box = box;
        numAtoms = box.getMoleculeList().size();
        dim = box.getSpace().D();
        rdot = new Vector[dim*nBeads*numAtoms];
        for (int i = 0; i < dim*nBeads*numAtoms; i++) {
            rdot[i] = box.getSpace().makeVector();
        }
        latticePositions = box.getSpace().makeVectorArray(box.getMoleculeList().size());
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

        getHOStagingParams();

        if (omegaSample != 0) {
            double En_ho_stage = hbar*omega/2.0*cothR;
            double Cvn_ho_stage = dim * nBeads / 2.0 * temperature * temperature;
            for (int k = 0; k < nBeads; k++) {
                En_ho_stage -= dim * gamma[k];
                Cvn_ho_stage += dim * dGamma[k];
            }
            Cvn_ho_stage *= beta*beta;
            System.out.println(" En_ho_stage:  " + En_ho_stage);
            System.out.println(" Cvn_ho_stage: " + Cvn_ho_stage);
        }
        massRing = nBeads*box.getLeafList().get(0).getType().getMass();

    }

    public void setNumShifts(int nShifts) {
        if (nShifts < 0) throw new RuntimeException("nShifts cannot be negative");
        this.nShifts = nShifts;
    }

    @Override
    public IData getData() {
        double[] x = data.getData();
        IMoleculeList molecules = box.getMoleculeList();

        pmBonding.computeAll(true);
        pcP1harm.computeAll(true); //no Cv
        pcP1ah.computeAll(true); //no Cv
//        pcP1harm.computeAll(true, this); //with Cv (replace 'this' by 'null' for no Cv)
//        pcP1ah.computeAll(true, this); //with Cv (replace 'this' by 'null' for no Cv)

        int numAtoms = molecules.size();
        double En0 = dim*nBeads*numAtoms/2.0/beta*R*cothR + R/sinhR*pcP1harm.getLastEnergy() + pcP1ah.getLastEnergy() - R*cothR*pmBonding.getLastEnergy();
        if (box.getMoleculeList().size() > 1) {
            En0 -= dim/2.0/beta;
        }
        double Cvn0 = nBeads/2.0/beta/beta - 2.0*pmBonding.getLastEnergy()/beta;
        for (int i=0; i<nBeads; i++) {
            En0 -= dim*numAtoms*gamma[i];
            Cvn0 += dim*numAtoms*dGamma[i];
        }

        Vector[] forcesUharm = pcP1harm.getForces();
        Vector[] forcesUah = pcP1ah.getForces();
        Vector[] forcesK = pmBonding.getForces();
        Vector tmp_r = box.getSpace().makeVector();
        Vector drj = box.getSpace().makeVector();
        Vector v = box.getSpace().makeVector();
        Vector a = box.getSpace().makeVector();
        Vector v0 = box.getSpace().makeVector();
        Vector a0 = box.getSpace().makeVector();
        Vector vPrev = box.getSpace().makeVector();
        Vector aPrev = box.getSpace().makeVector();

        Vector drShift = box.getSpace().makeVector();
        if (box.getMoleculeList().size() > 1) {
            drShift = computeShift();
        }
        int ns = nShifts+1;

        for (int i=0; i < x.length; i++){
            x[i] = 0;
        }

        for (int i = 0; i < numAtoms; i++) {
            IAtomList beads = molecules.get(i).getChildList();
            if (ns > beads.size() || (beads.size() / ns) * ns != beads.size()) {
                throw new RuntimeException("# of beads must be a multiple of (# of shifts + 1)");
            }
            for (int indexShift=0; indexShift<beads.size(); indexShift += beads.size()/ns) {
                double En = 0;
                double Cvn = 0;
                Vector dr0 = box.getSpace().makeVector();
                dr0.Ev1Mv2(beads.get(indexShift).getPosition(), latticePositions[i]);
                dr0.PE(drShift);
                box.getBoundary().nearestImage(dr0);

                Vector drPrev = box.getSpace().makeVector();

                for (int j = 0; j < beads.size(); j++) {
                    int aj = (j + indexShift) % beads.size();
                    IAtom atomj = beads.get(aj);
                    drj.Ev1Mv2(atomj.getPosition(), latticePositions[i]);
                    drj.PE(drShift);
                    box.getBoundary().nearestImage(drj);
                    tmp_r.E(drj);
                    if (j > 0) {
                        tmp_r.PEa1Tv1(-f11[j], drPrev);
                        tmp_r.PEa1Tv1(-f1N[j], dr0);
                    }
                    v.Ea1Tv1(gamma[j], tmp_r);
                    a.Ea1Tv1(dGamma[j]+gamma[j]*gamma[j], tmp_r);
                    if (j > 0) {
                        v.PEa1Tv1(df11[j], drPrev);
                        v.PEa1Tv1(df1N[j], dr0);
                        v.PEa1Tv1(f11[j], vPrev);
                        v.PEa1Tv1(f1N[j], v0);

                        a.PEa1Tv1(d2f11[j], drPrev);
                        a.PEa1Tv1(d2f1N[j], dr0);
                        a.PEa1Tv1(2.0*df11[j], vPrev);
                        a.PEa1Tv1(2.0*df1N[j], v0);
                        a.PEa1Tv1(f11[j], aPrev);
                        a.PEa1Tv1(f1N[j], a0);

                        Vector tmpV = box.getSpace().makeVector();
                        tmpV.Ev1Mv2(v, vPrev);
//                        Cvn -= beta*massRing*omegaN*omegaN*(tmpV.squared());
//                        if (j == nBeads-1){
//                            tmpV.Ev1Mv2(v0, v);
//                            Cvn -= beta*massRing*omegaN*omegaN*(tmpV.squared());
//                        }
                    }
                    int jj = atomj.getLeafIndex();
                    En -= beta*(forcesUharm[jj].dot(v) + forcesUah[jj].dot(v));
                    Cvn += beta*(forcesUharm[jj].dot(a)+forcesUah[jj].dot(a)) + 2.0*(forcesUharm[jj].dot(v)+forcesUah[jj].dot(v));
                    if (nBeads > 1) {
                        En -= beta*forcesK[jj].dot(v);
                        Cvn += beta*forcesK[jj].dot(a) - 2.0*forcesK[jj].dot(v);
                    }

                    drPrev.E(drj);
                    if (j == 0) {
                        v0.E(v);
                        a0.E(a);
                    }
                    vPrev.E(v);
                    aPrev.E(a);
                    rdot[jj].E(v);
                } // beads
                drdotHdrdot = 0;
//                pcP1.computeAll(false, this); // compute Hessian, using just-computed rdot
                Cvn -= beta*drdotHdrdot;
                x[0] += En0 + En;
                x[1] += Cvn0 + Cvn + (En0+En)*(En0+En);
            }//shifts
        }//atoms
        x[0] /= ns;
        x[1] /= ns;
        return data;
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector tmpV = box.getSpace().makeVector();
        tmpV.E(rdot[j]);
        Hij.transform(tmpV);
        drdotHdrdot += rdot[i].dot(tmpV); //drdot.H.drdot
    }

    public boolean wantsHessian() {
        return true;
    }

    protected Vector computeShift() {
        if (box.getMoleculeList().size() == 1) {
            return box.getSpace().makeVector();
        }
        int n = box.getMoleculeList().get(0).getChildList().size();
        Vector shift0 = box.getSpace().makeVector();
        Boundary boundary = box.getBoundary();
        Vector dr = box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        dr.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[0]);
        boundary.nearestImage(dr);
        shift0.Ea1Tv1(-1, dr);
        // will shift ring0 back to lattice site; everything should be close and PBC should lock in
        // now determine additional shift needed to bring back to original COM
        Vector totalShift = box.getSpace().makeVector();
        for (int j = 0; j < box.getMoleculeList().size(); j++) {
            IMolecule m = box.getMoleculeList().get(j);
            for (int i = 0; i < n; i++) {
                Vector r = m.getChildList().get(i).getPosition();
                dr.Ev1Mv2(r, latticePositions[j]);
                dr.PE(shift0);
                boundary.nearestImage(dr);
                totalShift.PE(dr);
            }
        }
        totalShift.TE(-1.0/box.getLeafList().size());
        totalShift.PE(shift0);
        return totalShift;
    }



    public  void getHOStagingParams(){
        R = beta*hbar*omega/nBeads;
        double dR = hbar*omega/nBeads;
        double dR2 = dR*dR;
        double d2R = 0;
        sinhR = Math.sinh(R);
        double coshR = Math.cosh(R);
        cothR = coshR/sinhR;
        double sinhNR = Math.sinh(nBeads*R);
        double coshhNR = Math.cosh(nBeads*R);

        f11 = new double[nBeads];
        f1N = new double[nBeads];
        df11 = new double[nBeads];
        df1N = new double[nBeads];
        d2f11 = new double[nBeads];
        d2f1N = new double[nBeads];
        gamma = new double[nBeads];
        dGamma = new double[nBeads];

        gamma[0] = omegaSample == 0 ? 0 : -hbar*omega/2/Math.sinh(beta*hbar*omega);
        dGamma[0] = omegaSample == 0 ? 0 : 111111;

        for (int i=1; i<nBeads; i++) {
            double sinhNmiA = Math.sinh((nBeads-i)*R);
            double coshNmiA = Math.cosh((nBeads-i)*R);
            double sinhNmip1A = Math.sinh((nBeads-i+1)*R);
            double coshNmip1A = Math.cosh((nBeads-i+1)*R);

            gamma[i] = omegaSample == 0 ? 1.0/2.0/beta : -hbar*omega/2.0/nBeads*((nBeads-i+1)/Math.tanh((nBeads-i+1)*R)
                    -(nBeads-i)/Math.tanh((nBeads-i)*R) - 1.0/Math.tanh(R));
            dGamma[i] = omegaSample == 0 ? -1.0/2.0/beta/beta : 111111;
            f11[i] = omegaSample == 0 ? (nBeads-i)/(nBeads-i+1.0) : sinhNmiA/sinhNmip1A;
            f1N[i] = omegaSample == 0 ? 1.0/(nBeads-i+1) : sinhR/sinhNmip1A;

            df11[i] = omegaSample == 0 ? 0 : dR/sinhNmip1A*((nBeads-i)*coshNmiA-(nBeads-i+1)*sinhNmiA*coshNmip1A/sinhNmip1A);
            df1N[i] = omegaSample == 0 ? 0 : dR/sinhNmip1A*(coshR - (nBeads-i+1)*sinhR*coshNmip1A/sinhNmip1A);

            d2f11[i] = omegaSample == 0 ? 0 : d2R/dR*df11[i] + dR2/sinhNmip1A*((nBeads-i)*(nBeads-i)*sinhNmiA
                    - 2*(nBeads-i+1)*(nBeads-i)*coshNmiA*coshNmip1A/sinhNmip1A
                    + 0.5*(nBeads-i+1)*(nBeads-i+1)*sinhNmiA/sinhNmip1A/sinhNmip1A*(3.0+Math.cosh(2*(nBeads-i+1)*R)));
            d2f1N[i] = omegaSample == 0 ? 0 : d2R/dR*df1N[i] + dR2/sinhNmip1A*(sinhR
                    -2*(nBeads-i+1)*coshR*coshNmip1A/sinhNmip1A
                    + 0.5*(nBeads-i+1)*(nBeads-i+1)*sinhR/sinhNmip1A/sinhNmip1A*(3.0 + Math.cosh(2*(nBeads-i+1)*R)));
        }

    }
    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }
}
