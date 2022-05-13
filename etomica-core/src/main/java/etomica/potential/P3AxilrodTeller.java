/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.simulation.Simulation;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.Space3D;
import etomica.species.ISpecies;
import etomica.species.SpeciesGeneral;
import etomica.util.random.RandomMersenneTwister;

import java.util.HashMap;
import java.util.Map;

/**
 * Axilrod-Teller potential.  The potential is atomic.  Ionization energy and
 * polarizability are required as input for each atom type. 
 * 
 * @author Andrew Schultz
 */
public class P3AxilrodTeller implements IPotential3 {

    protected final Map<AtomType, MyAgent> paramsManager;
    protected double[][][] epsilon;
    protected final Space space;
    protected final Vector dr1, dr2;
    protected Boundary boundary;
    protected double range;
    
    public P3AxilrodTeller(Space space, Map<AtomType, MyAgent> paramsManager, double rc) {
        this.space = space;
        this.paramsManager = paramsManager;
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        range = rc;
        epsilon = new double[0][0][0];
    }

    private void copyEpsilonUp(int m) {
        if (epsilon.length <= m) {
            double[][][] e = new double[m+1][m+1][m+1];
            for (int i=0; i<epsilon.length; i++) {
                for (int j=0; j<epsilon.length; j++) {
                    for (int k=0; k<epsilon.length; k++) {
                        e[i][j][k] = epsilon[i][j][k];
                    }
                }
            }
            epsilon = e;
        }
    }

    public void setEpsilon(AtomType type1, AtomType type2, AtomType type3, double eps) {
        int m = type3.getIndex();
        if (epsilon.length <= m) {
            copyEpsilonUp(m);
        }
        double[][][] e = epsilon;
        for (int i=0; i<epsilon.length; i++) {
            for (int j = 0; j < epsilon.length; j++) {
                for (int k = 0; k < epsilon.length; k++) {
                    e[k][j][i] = e[j][k][i] = e[j][i][k] = e[k][i][j] = e[i][k][j] = e[i][j][k] = eps;
                }
            }
        }
    }

    @Override
    public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3, double[] virial) {

        double RAB2 = dr12.squared();
        double RAC2 = dr13.squared();
        double RBC2 = dr23.squared();
        if (RAB2*RAC2*RBC2 == 0) return 0;

        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);
        if (RAB>range || RAC>range || RBC>range) return 0;

        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);
        double cp = 3 * costhetaA * costhetaB * costhetaC;

        int m = Math.max(Math.max(atom1.getType().getIndex(), atom2.getType().getIndex()), atom3.getType().getIndex());
        if (epsilon.length <= m) {
            copyEpsilonUp(m);
        }

        if (epsilon[atom1.getType().getIndex()][atom2.getType().getIndex()][atom3.getType().getIndex()] != 0) {
            double eps = epsilon[atom1.getType().getIndex()][atom2.getType().getIndex()][atom3.getType().getIndex()];
            double u = eps*(cp+1)/(RAB*RAB2*RAC*RAC2*RBC*RBC2);
            virial[0] = -9*u;
            return u;
        }
        IAtom[] atoms = new IAtom[]{atom1, atom2, atom3};
        double ep = 1;
        double es = 0;
        double[] eps = new double[3];
        double ap = 1;
        for (int i=0; i<3; i++) {
            MyAgent ag = paramsManager.get(atoms[i].getType());
            eps[(i+1)%3] += ag.E;
            eps[(i+2)%3] += ag.E;
            ep *= ag.E;
            es += ag.E;

            ap *= ag.alpha;
        }
        double e123 = ep*es;
        for (int i=0; i<3; i++) {
            e123 /= eps[i];
        }

        double u = 1.5*e123*ap*(cp+1)/(RAB*RAB2*RAC*RAC2*RBC*RBC2);
        virial[0] = -9*u;
        return u;
    }

    @Override
    public double udu(Vector drAB, Vector drAC, Vector drBC, IAtom atomA, IAtom atomB, IAtom atomC, double[] virial, Vector fA, Vector fB, Vector fC) {
        double RAB2 = drAB.squared();
        double RAC2 = drAC.squared();
        double RBC2 = drBC.squared();
        if (RAB2*RAC2*RBC2 == 0) return 0;

        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);
        if (RAB>range || RAC>range || RBC>range) return 0;

        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        if (costhetaA == 0) costhetaA = 1e-9;
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        if (costhetaB == 0) costhetaB = 1e-9;
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);
        if (costhetaC == 0) costhetaC = 1e-9;
        double cp = 3 * costhetaA * costhetaB * costhetaC;

        double myEpsilon = 0;
        double epsCheck = epsilon[atomA.getType().getIndex()][atomB.getType().getIndex()][atomB.getType().getIndex()];
        if (epsCheck != 0) {
            myEpsilon = epsCheck;
        }
        else {
            IAtom[] atoms = new IAtom[]{atomA, atomB, atomC};
            double ep = 1;
            double es = 0;
            double[] eps = new double[3];
            double ap = 1;
            for (int i = 0; i < 3; i++) {
                MyAgent ag = paramsManager.get(atoms[i].getType());
                eps[(i + 1) % 3] += ag.E;
                eps[(i + 2) % 3] += ag.E;
                ep *= ag.E;
                es += ag.E;

                ap *= ag.alpha;
            }
            double e123 = ep * es;
            for (int i = 0; i < 3; i++) {
                e123 /= eps[i];
            }
            myEpsilon = 1.5*e123*ap;
        }

        double u = myEpsilon*(cp+1)/(RAB*RAB2*RAC*RAC2*RBC*RBC2);

        double rdudr = -3*u;
        double cosdudcos = myEpsilon*cp/(RAB*RAB2*RAC*RAC2*RBC*RBC2);

        // fA contribution
        Vector tmp = Vector.d(drAB.getD());
        tmp.Ea1Tv1(-rdudr/RAB2, drAB);
        tmp.PEa1Tv1(-rdudr/RAC2, drAC);

        Vector tmp2 = Vector.d(drAB.getD());  // dcosthetaBdrA
        tmp2.Ea1Tv1(1/RAB/RBC, drBC);
        tmp2.PEa1Tv1(costhetaB/RAB2, drAB);
        Vector dcosthetaBdrA = Vector.d(drAB.getD());
        dcosthetaBdrA.E(tmp2);

        tmp.PEa1Tv1(cosdudcos/costhetaB, tmp2);

        // dcosthetaCdrA
        tmp2.Ea1Tv1(1/RAC/RBC, drBC);
        tmp2.PEa1Tv1(-costhetaC/RAC2, drAC);
        tmp.PEa1Tv1(-cosdudcos/costhetaC, tmp2);

        // dcothetaAdrA = -(dcosthetaAdrB + dcosthetaAdrC)
        tmp2.Ea1Tv1(1/RAB/RAC, drAC); // dcosthetaAdrB
        tmp2.PEa1Tv1(-costhetaA/RAB2, drAB);
        Vector dcosthetaAdrB = Vector.d(drAB.getD());
        dcosthetaAdrB.E(tmp2);
        tmp2.PEa1Tv1(1/RAC/RAB, drAB); // dcosthetaAdrC
        tmp2.PEa1Tv1(-costhetaA/RAC2, drAC);
        tmp.PEa1Tv1(-cosdudcos/costhetaA, tmp2);
        // our tmp is the gradient for A, so subtract to get force
        fA.ME(tmp);
        fC.PE(tmp);

        // fB
        tmp.Ea1Tv1(rdudr/RAB2, drAB);
        tmp.PEa1Tv1(-rdudr/RBC2, drBC);

        // dcosthetaAdrB
        tmp.PEa1Tv1(cosdudcos/costhetaA, dcosthetaAdrB);

        // dcosthetaCdrB
        tmp2.Ea1Tv1(1/RBC/RAC, drAC);
        tmp2.PEa1Tv1(-costhetaC/RBC2, drBC);
        tmp.PEa1Tv1(-cosdudcos/costhetaC, tmp2);

        // dcothetaBdrB = -(dcosthetaBdrA + dcosthetaBdrC)
        tmp2.E(dcosthetaBdrA); // dcosthetaBdrA
        tmp2.PEa1Tv1(-1/RBC/RAB, drAB); // dcosthetaBdrC
        tmp2.PEa1Tv1(-costhetaB/RBC2, drBC);
        tmp.PEa1Tv1(-cosdudcos/costhetaB, tmp2);
        // our tmp is the gradient for B, so subtract to get force
        fB.ME(tmp);
        fC.PE(tmp);

        virial[0] = -9*u;
        return u;

    }

    public double getRange() {
        return range;
    }

    // CO2: alpha=2.913, alpha(omega)=2.639, I=13.7eV
    // H2O: alpha=1.444, I=12.6eV
    public static class MyAgent {
        public final double alpha, E;
        /**
         * @params alpha polarizability
         * @params E ionization energy
         */
        public MyAgent(double alpha, double E) {
            this.alpha = alpha;
            this.E = E;
        }
    }

    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        sim.setRandom(new RandomMersenneTwister(2));
        Space space = sim.getSpace();
        ISpecies species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(sim));
        sim.addSpecies(species);
        Box box = sim.makeBox();
        box.setNMolecules(species, 3);
        IAtomList atoms = box.getLeafList();
        Map<AtomType,MyAgent> agentMap = new HashMap<>();
        agentMap.put(species.getLeafType(), new MyAgent(1, 1));
        P3AxilrodTeller p3 = new P3AxilrodTeller(space, agentMap, Double.POSITIVE_INFINITY);
        Vector r0 = atoms.get(0).getPosition();
        Vector r1 = atoms.get(1).getPosition();
        r1.setRandomInSphere(sim.getRandom());
        r1.TE(2);
        Vector r2 = atoms.get(2).getPosition();
        r2.setRandomInSphere(sim.getRandom());
        r2.TE(2);
        Vector drAB = space.makeVector();
        Vector drAC = space.makeVector();
        Vector drBC = space.makeVector();
        drAB.E(r1);
        drAC.E(r2);
        drBC.Ev1Mv2(r2, r1);
        Vector fA = space.makeVector(), fB = space.makeVector(), fC = space.makeVector();
        double u0 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);

        double del = 0.00001;
        Vector dr = space.makeVector();
        dr.setRandomSphere(sim.getRandom());
        dr.TE(del);
        r2.PE(dr);
        drAB.Ev1Mv2(r1, r0);
        drAC.Ev1Mv2(r2, r0);
        drBC.Ev1Mv2(r2, r1);
        double u1 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);
        fA.TE(0.5);
        fB.TE(0.5);
        fC.TE(0.5);

        double fdr = fC.dot(dr);
        System.out.println(-fdr+" "+(u1-u0));

    }
}
