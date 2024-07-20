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

public class MeterPIHMAReal2 implements IDataSource, PotentialCallback {
    protected final PotentialMasterBonding pmBonding;
    protected final PotentialCompute pcP1;
    protected double beta;
    protected final DataTag tag;
    protected DataDoubleArray.DataInfoDoubleArray dataInfo;
    protected DataDoubleArray data;
    protected final MCMoveHOReal2 move;
    protected int nShifts = 0;
    protected int dim;
    protected Vector[] rdot;
    protected double drdotHdrdot;
    protected Vector[] latticePositions;
    protected int numAtoms, nBeads;
    protected double[] gamma, dGamma, f11, f1N, df11, df1N, d2f11, d2f1N;
    protected double[][] centerCoefficients, dCenterCoefficients, d2CenterCoefficients;
    protected Box box;
    protected double EnShift;
    protected double omegaN, mass;
    protected Vector drShift;


    public MeterPIHMAReal2(PotentialMasterBonding pmBonding, PotentialCompute pcP1, int nBeads, double temperature, MCMoveHOReal2 move, int nShifts) {
        int nData = 2;
        this.move = move;
        this.beta = 1/temperature;
        this.nBeads = nBeads;
        data = new DataDoubleArray(nData);
        dataInfo = new DataDoubleArray.DataInfoDoubleArray("PI",Null.DIMENSION, new int[]{nData});
        tag = new DataTag();
        dataInfo.addTag(tag);
        this.pmBonding = pmBonding;
        this.pcP1 = pcP1;
        this.box = move.getBox();
        this.omegaN = nBeads/(beta*move.hbar);
        numAtoms = box.getMoleculeList().size();
        this.mass = box.getMoleculeList().get(0).getType().getMass(); //RP mass
        dim = box.getSpace().D();
        rdot = box.getSpace().makeVectorArray(numAtoms*nBeads);
        latticePositions = box.getSpace().makeVectorArray(numAtoms);
        for (int i=0; i<latticePositions.length; i++) {
            latticePositions[i].E(CenterOfMass.position(box, box.getMoleculeList().get(i)));
        }

        gamma = move.getGamma();
        dGamma = move.getdGamma();
        centerCoefficients = move.getCenterCoefficients();
        f11 = centerCoefficients[0];
        f1N = centerCoefficients[1];
        dCenterCoefficients = move.getDCenterCoefficients();
        df11 = dCenterCoefficients[0];
        df1N = dCenterCoefficients[1];
        d2CenterCoefficients = move.getD2CenterCoefficients();
        d2f11 = d2CenterCoefficients[0];
        d2f1N = d2CenterCoefficients[1];

        if (move.omega2 != 0) {
            double En_ho_stage = dim*nBeads/2.0/beta;
            double Cvn_ho_stage = dim*nBeads/2.0*temperature * temperature;
            for (int k = 0; k < nBeads; k++) {
                En_ho_stage -= dim * gamma[k];
                Cvn_ho_stage += dim * dGamma[k];
            }
            Cvn_ho_stage *= beta * beta;
            //COM
            En_ho_stage -= dim/2.0/beta;
            Cvn_ho_stage -= dim/2.0/beta/beta;
            System.out.println(" En_ho_stage:  " + En_ho_stage);
            System.out.println(" Cvn_ho_stage: " + Cvn_ho_stage);
        }

        this.EnShift = 0;
        this.nShifts = nShifts;
        drShift = box.getSpace().makeVector();

        if (this.nShifts < 0) throw new RuntimeException("nShifts cannot be negative");
        int ns = this.nShifts + 1;
        if (ns > nBeads || (nBeads / ns) * ns != nBeads) {
            throw new RuntimeException("# of beads must be a multiple of (# of shifts + 1)");
        }
    }

    @Override
    public IData getData() {
        Box box = move.getBox();
        double[] x = data.getData();
        for (int i = 0; i < x.length; i++) x[i] = 0;

        pmBonding.computeAll(true); // only forces
        pcP1.computeAll(true); // only forces

        IMoleculeList molecules = box.getMoleculeList();
        double En0 = dim*numAtoms*nBeads/2.0/beta + pcP1.getLastEnergy() - pmBonding.getLastEnergy();
        double Cvn0 = dim*numAtoms*nBeads/2.0/beta/beta - 2.0/beta*pmBonding.getLastEnergy();
        for (int i=0; i<nBeads; i++) {
            En0 -= dim*numAtoms*gamma[i];
            Cvn0 += dim*numAtoms*dGamma[i];
        }

        if (numAtoms > 1 && move.omega2 != 0) {
            En0 -= dim/2.0/beta;
            Cvn0 -= dim/2.0/beta/beta;
        }

        Vector[] forcesU = pcP1.getForces();
        Vector[] forcesK = pmBonding.getForces();
        Vector tmp_r = box.getSpace().makeVector();
        Vector drj = box.getSpace().makeVector();
        Vector v = box.getSpace().makeVector();
        Vector a = box.getSpace().makeVector();
        Vector v0 = box.getSpace().makeVector();
        Vector a0 = box.getSpace().makeVector();
        Vector vPrev = box.getSpace().makeVector();
        Vector aPrev = box.getSpace().makeVector();

        if (box.getMoleculeList().size() > 1) drShift = computeShift();
        int ns = nShifts+1;
        Vector dr0 = box.getSpace().makeVector();


        double y = 0;
        double U = pcP1.getLastEnergy();

        for (int indexShift = 0; indexShift < nBeads; indexShift += nBeads/ns) {
            double En = 0;
            double Cvn = 0;
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
                        Cvn -= beta*mass/nBeads*omegaN*omegaN*(tmpV.squared()); //mass*move is RP mass (=1)
                        if (j == nBeads-1){
                            tmpV.Ev1Mv2(v0, v);
                            Cvn -= beta*mass/nBeads*omegaN*omegaN*(tmpV.squared());
                        }
                    }
                    int jj = atomj.getLeafIndex();
//                    En -= beta*forcesU[jj].dot(v) + beta*forcesK[jj].dot(v);
                    y -= beta*forcesU[jj].dot(v);
                    Cvn += beta*(forcesU[jj].dot(a) + forcesK[jj].dot(a)) + 2.0*(forcesU[jj].dot(v) - forcesK[jj].dot(v));

                    drPrev.E(drj);
                    if (j == 0) {
                        v0.E(v);
                        a0.E(a);
                    }
                    vPrev.E(v);
                    aPrev.E(a);
                    rdot[jj].E(v);
                } // beads
            }//mol
            drdotHdrdot = 0;
            pcP1.computeAll(false, this); // compute Hessian, using just-computed rdot
            Cvn -= beta*drdotHdrdot;

            x[0] += En;
            x[1] += Cvn + (En0+En- EnShift)*(En0+En- EnShift);

//            double y =  -2.0/beta*pmBonding.getLastEnergy() + Cvn;
        }//shifts
        x[0] = En0 + x[0]/ns - EnShift;
        x[1] = Cvn0 + x[1]/ns;

        y /= ns;

//        y += pcP1.getLastEnergy() - pmBonding.getLastEnergy();
        y += U;

        System.out.println("AN: " + y);

        return data;
    }

    public void fieldComputeHessian(int i, Tensor Hii) {
        Vector tmpV = move.getBox().getSpace().makeVector();
        tmpV.E(rdot[i]);
        Hii.transform(tmpV);
        drdotHdrdot += rdot[i].dot(tmpV); //drdot.H.drdot
    }

    public void pairComputeHessian(int i, int j, Tensor Hij) { // in general potential, Hij is the Hessian between same beads of atom i and j
        Vector rdotij = box.getSpace().makeVector();
        rdotij.Ev1Mv2(rdot[i], rdot[j]);
        Hij.transform(rdotij);
        drdotHdrdot += rdotij.dot(rdot[i]) - rdotij.dot(rdot[j]);
    }

    public boolean wantsHessian() {
        return true;
    }

    public IDataInfo getDataInfo() {
        return dataInfo;
    }

    public DataTag getTag() {
        return tag;
    }

    public void setEnShift(double E) { this.EnShift = E; }

    protected Vector computeShift() {
        if (box.getMoleculeList().size() == 1) return box.getSpace().makeVector();
        IAtomList atoms = box.getMoleculeList().get(0).getChildList();
        Vector drShift = box.getSpace().makeVector();
        drShift.Ev1Mv2(atoms.get(0).getPosition(), latticePositions[0]);
        drShift.TE(-1);
        return drShift;
    }
}