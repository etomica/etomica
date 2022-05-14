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
import etomica.species.SpeciesManager;
import etomica.units.Erg;

/**
 * 3-body contribution to Stillinger-Weber potential.
 *
 * <a href="https://doi.org/10.1103/PhysRevB.31.5262">DOI:10.1103/PhysRevB.31.5262</a>
 */
public class P3StillingerWeber implements IPotential3 {

    protected final double[][][] lambda, costheta0, epsilon;
    protected final double[][] gamma, a, sigma;
    protected final Space space;
    protected final Vector dr1, dr2;
    protected Boundary boundary;
    protected double range, range2;

    public P3StillingerWeber(SpeciesManager sm, Space space, double rc) {
        this.space = space;
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        range = rc;
        range2 = range*range;
        int nat = sm.getAtomTypeCount();
        lambda = new double[nat][nat][nat];
        costheta0 = new double [nat][nat][nat];
        epsilon = new double[nat][nat][nat];
        gamma = new double[nat][nat];
        a = new double[nat][nat];
        sigma = new double[nat][nat];
    }

    public void setParameters2(AtomType type1, AtomType type2, double gamma, double a, double sigma) {
        this.gamma[type1.getIndex()][type2.getIndex()] = this.gamma[type2.getIndex()][type1.getIndex()] = gamma;
        this.a[type1.getIndex()][type2.getIndex()] = this.a[type2.getIndex()][type1.getIndex()] = a;
        this.sigma[type1.getIndex()][type2.getIndex()] = this.sigma[type2.getIndex()][type1.getIndex()] = sigma;
    }

    public void setParameters3(AtomType type1, AtomType type2, AtomType type3, double lambda, double costheta0, double epsilon) {
        this.lambda[type1.getIndex()][type2.getIndex()][type3.getIndex()] =
                this.lambda[type3.getIndex()][type2.getIndex()][type1.getIndex()] = lambda;
        this.epsilon[type1.getIndex()][type2.getIndex()][type3.getIndex()] =
                this.lambda[type3.getIndex()][type2.getIndex()][type1.getIndex()] = epsilon;
        this.costheta0[type1.getIndex()][type2.getIndex()][type3.getIndex()] =
                this.costheta0[type3.getIndex()][type2.getIndex()][type1.getIndex()] = costheta0;
    }

    @Override
    public double u(Vector dr12, Vector dr13, Vector dr23, IAtom atom1, IAtom atom2, IAtom atom3, double[] virial) {

        // we get called once with this triplet.  1(A) is our central atom
        double RAB2 = dr12.squared();
        double RAC2 = dr13.squared();
        double RBC2 = dr23.squared();
        if (RAB2*RAC2*RBC2 == 0 || RAB2 > range2 || RAC2 > range2) return 0;

        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);

        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);

        int t1 = atom1.getType().getIndex();
        int t2 = atom2.getType().getIndex();
        int t3 = atom3.getType().getIndex();
        double l1 = lambda[t2][t1][t3];
        double e1 = epsilon[t2][t1][t3];
        double c01 = costheta0[t2][t1][t3];
        double g12 = gamma[t1][t2], g13 = gamma[t1][t3];
        double a12 = a[t1][t2], a13 = a[t1][t3];
        double s12 = sigma[t1][t2], s13 = sigma[t1][t3];
        double arg = g12/(RAB/s12-a12);
        double exp12 = Math.exp(arg);
        double r12dexp12dr12 = -exp12*arg*RAB/(RAB/s12-a12);
        arg = g13/(RAC/s13-a13);
        double exp13 = Math.exp(arg);
        double r13dexp13dr13 = -exp13*arg*RAC/(RAC/s13-a13);
        double deltacos = costhetaA - c01;
        double u1 = e1*l1*deltacos*deltacos*exp12*exp13;

        double v1 = u1*(1/exp12*r12dexp12dr12 + 1/exp13*r13dexp13dr13);

        virial[0] = v1;
        return u1;
    }

    @Override
    public double udu(Vector drAB, Vector drAC, Vector drBC, IAtom atomA, IAtom atomB, IAtom atomC, double[] virial, Vector fA, Vector fB, Vector fC) {
        double RAB2 = drAB.squared();
        double RAC2 = drAC.squared();
        double RBC2 = drBC.squared();
        if (RAB2*RAC2*RBC2 == 0 || RAB2 > range2 || RAC2 > range2) return 0;

        double RAB = Math.sqrt(RAB2);
        double RAC = Math.sqrt(RAC2);
        double RBC = Math.sqrt(RBC2);

        double costhetaA = (RAB2 + RAC2 - RBC2)/(2*RAC*RAB);
        if (costhetaA == 0) costhetaA = 1e-9;
        // costhetaB and C aren't going to be used by forceHelper but they're easy enough to compute here
        double costhetaB = (RAB2 + RBC2 - RAC2)/(2*RAB*RBC);
        if (costhetaB == 0) costhetaB = 1e-9;
        double costhetaC = (RAC2 + RBC2 - RAB2)/(2*RAC*RBC);
        if (costhetaC == 0) costhetaC = 1e-9;

        int t1 = atomA.getType().getIndex();
        int t2 = atomB.getType().getIndex();
        int t3 = atomC.getType().getIndex();
        double l1 = lambda[t2][t1][t3];
        double e1 = epsilon[t2][t1][t3];
        double c01 = costheta0[t2][t1][t3];
        double g12 = gamma[t1][t2], g13 = gamma[t1][t3];
        double a12 = a[t1][t2], a13 = a[t1][t3];
        double s12 = sigma[t1][t2], s13 = sigma[t1][t3];
        double arg = g12/(RAB/s12-a12);
        double exp12 = arg < 0 ? Math.exp(arg) : 0;
        double r12dexp12dr12 = -exp12*arg*RAB/s12/(RAB/s12-a12);
        arg = g13/(RAC/s13-a13);
        double exp13 = arg < 0 ? Math.exp(arg) : 0;
        double r13dexp13dr13 = -exp13*arg*RAC/s13/(RAC/s13-a13);

        double deltacos = costhetaA - c01;
        double u1 = e1*l1*deltacos*deltacos*exp12*exp13;
        double cosdu1dcos = deltacos!=0 ? 2*u1/deltacos*costhetaA : 0;

        double r12du1dr12 = exp12>0 ? u1/exp12*r12dexp12dr12 : 0, r13du1dr13 = exp13>0 ? u1/exp13*r13dexp13dr13 : 0;
        double v1 = r12du1dr12 + r13du1dr13;

        forceHelper(drAB, drAC, drBC, RAB, RAC, RBC, costhetaA, costhetaB, costhetaC, r12du1dr12, r13du1dr13,0, cosdu1dcos,0, 0, fA, fB, fC);
        if (fA.isNaN() || fB.isNaN() || fC.isNaN()) {
            throw new RuntimeException("oops "+fA+" "+fB+" "+fC);
        }

        virial[0] = v1;
        return u1;
    }

    public double getRange() {
        return range;
    }

    public static void main(String[] args) {
        Simulation sim = new Simulation(Space3D.getInstance());
        Space space = sim.getSpace();
        ISpecies species = SpeciesGeneral.monatomic(space, AtomType.simpleFromSim(sim));
        sim.addSpecies(species);
        Box box = sim.makeBox();
        box.setNMolecules(species, 3);
        IAtomList atoms = box.getLeafList();
        P3StillingerWeber p3 = new P3StillingerWeber(sim.getSpeciesManager(), sim.getSpace(), 1.8*2.0951);
        AtomType leafType = species.getLeafType();
        p3.setParameters2(leafType, leafType, 1.2, 1.8, 2.0951);
        p3.setParameters3(leafType, leafType, leafType, 21, -1.0/3.0, Erg.UNIT.toSim(3.4723e-12));
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
        double u00 = p3.u(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1]);
        double u0 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);
        System.out.println(u00+" "+u0);

        double del = 0.00001;
        Vector dr = space.makeVector();
        dr.setRandomSphere(sim.getRandom());
        dr.TE(del);
        r0.PE(dr);
        drAB.Ev1Mv2(r1, r0);
        drAC.Ev1Mv2(r2, r0);
        drBC.Ev1Mv2(r2, r1);
        double u1 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);
        fA.TE(0.5);
        fB.TE(0.5);
        fC.TE(0.5);

        double fdr = fA.dot(dr);
        System.out.println(-fdr+" "+(u1-u0));

        fA.E(0);
        fB.E(0);
        fC.E(0);

        u0 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);

        dr.setRandomSphere(sim.getRandom());
        dr.TE(del);
        r1.PE(dr);
        drAB.Ev1Mv2(r1, r0);
        drAC.Ev1Mv2(r2, r0);
        drBC.Ev1Mv2(r2, r1);
        u1 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);
        fA.TE(0.5);
        fB.TE(0.5);
        fC.TE(0.5);

        fdr = fB.dot(dr);
        System.out.println(-fdr+" "+(u1-u0));

        fA.E(0);
        fB.E(0);
        fC.E(0);

        u0 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);

        dr.setRandomSphere(sim.getRandom());
        dr.TE(del);
        r2.PE(dr);
        drAB.Ev1Mv2(r1, r0);
        drAC.Ev1Mv2(r2, r0);
        drBC.Ev1Mv2(r2, r1);
        u1 = p3.udu(drAB, drAC, drBC, atoms.get(0), atoms.get(1), atoms.get(2), new double[1], fA, fB, fC);
        fA.TE(0.5);
        fB.TE(0.5);
        fC.TE(0.5);

        fdr = fC.dot(dr);
        System.out.println(-fdr+" "+(u1-u0));

    }
}
