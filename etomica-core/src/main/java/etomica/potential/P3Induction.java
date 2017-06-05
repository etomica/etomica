/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.space.Boundary;
import etomica.box.Box;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.IPotentialAtomic;
import etomica.space.Vector;
import etomica.atom.AtomPair;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtomOriented;
import etomica.atom.MoleculePair;
import etomica.chem.elements.ElementSimple;
import etomica.chem.elements.Hydrogen;
import etomica.chem.elements.Oxygen;
import etomica.models.water.ConformationWaterGCPM;
import etomica.models.water.P2WaterSzalewicz;
import etomica.models.water.P2WaterSzalewicz.Component;
import etomica.models.water.PNWaterGCPM;
import etomica.models.water.SpeciesWater4P;
import etomica.simulation.Simulation;
import etomica.space.IOrientation;
import etomica.space.Space;
import etomica.space3d.OrientationFull3D;
import etomica.space3d.Space3D;
import etomica.species.SpeciesSpheresRotating;
import etomica.util.random.RandomMersenneTwister;
import etomica.util.random.RandomNumberGeneratorUnix;

/**
 * 3-body induction potential based on form used by Oakley and Wheatley.
 * 
 * http://dx.doi.org/10.1063/1.3059008
 *
 * @author Andrew Schultz
 */
public class P3Induction implements IPotentialAtomic {

    protected final AtomTypeAgentManager paramsManager;
    protected final Space space;
    protected final double[] I = new double[3];
    protected final double[] alpha = new double[3];
    protected final Vector dr1, dr2;
    protected final Vector ri, rj, rk;
    protected final Vector rij, rik;
    protected final Vector or3;
    protected Boundary boundary;

    public P3Induction(Space space, AtomTypeAgentManager paramsManager) {
        this.space = space;
        this.paramsManager = paramsManager;
        dr1 = space.makeVector();
        dr2 = space.makeVector();
        ri = space.makeVector();
        rj = space.makeVector();
        rk = space.makeVector();
        rij = space.makeVector();
        rik = space.makeVector();
        or3 = space.makeVector();
    }

    public double energy(IAtomList atoms) {
        double sum = 0;
        for (int i=0; i<3; i++) {
            
            int j = (i+1)%3;
            int k = (i+2)%3;
            IAtomOriented atomi = (IAtomOriented)atoms.getAtom(i);
            IOrientation ori = atomi.getOrientation();
            IAtomOriented atomj = (IAtomOriented)atoms.getAtom(j);
            IOrientation orj = atomj.getOrientation();
            IAtomOriented atomk = (IAtomOriented)atoms.getAtom(k);
            IOrientation ork = atomk.getOrientation();
            MyAgent agi = (MyAgent)paramsManager.getAgent(atoms.getAtom(i).getType());
            MyAgent agj = (MyAgent)paramsManager.getAgent(atoms.getAtom(j).getType());
            MyAgent agk = (MyAgent)paramsManager.getAgent(atoms.getAtom(k).getType());
            for (int ip=0; ip<agi.alpha.length; ip++) {
                ri.E(atomi.getPosition());
                ri.PEa1Tv1(agi.polSite[ip].getX(0), ori.getDirection());
                if (ori instanceof OrientationFull3D) {
                    Vector or2 = ((OrientationFull3D)ori).getSecondaryDirection();
                    ri.PEa1Tv1(agi.polSite[ip].getX(1), or2);
                    or3.E(ori.getDirection());
                    or3.XE(or2);
                    ri.PEa1Tv1(agi.polSite[ip].getX(2), or3);
                }
                for (int jq=0; jq<agj.q.length; jq++) {
                    rj.E(atomj.getPosition());
                    rj.PEa1Tv1(agj.qSite[jq].getX(0), orj.getDirection());
                    if (orj instanceof OrientationFull3D) {
                        Vector or2 = ((OrientationFull3D)orj).getSecondaryDirection();
                        rj.PEa1Tv1(agj.qSite[jq].getX(1), or2);
                        or3.E(orj.getDirection());
                        or3.XE(or2);
                        rj.PEa1Tv1(agj.qSite[jq].getX(2), or3);
                    }
                    rij.Ev1Mv2(rj, ri);
                    double rij2 = rij.squared();
                    // normalize rij
                    rij.TE(1/Math.sqrt(rij2));
                    for (int kq=0; kq<agk.q.length; kq++) {
                        rk.E(atomk.getPosition());
                        rk.PEa1Tv1(agk.qSite[kq].getX(0), ork.getDirection());
                        if (ork instanceof OrientationFull3D) {
                            Vector or2 = ((OrientationFull3D)ork).getSecondaryDirection();
                            rk.PEa1Tv1(agk.qSite[kq].getX(1), or2);
                            or3.E(ork.getDirection());
                            or3.XE(or2);
                            rk.PEa1Tv1(agk.qSite[kq].getX(2), or3);
                        }
                        rik.Ev1Mv2(rk, ri);
                        double rik2 = rik.squared();
                        rik.TE(1/Math.sqrt(rik2));
                        
                        sum += agi.alpha[ip]*agj.q[jq]*agk.q[kq]/(rij2*rik2)*rij.dot(rik); 
                    }
                }
            }
        }

        return -sum;
    }

    public double getRange() {
        return Double.POSITIVE_INFINITY;
    }

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }

    public int nBody() {
        return 3;
    }

    public static class MyAgent {
        public final double[] alpha, q;
        public final Vector[] polSite, qSite;
        public MyAgent(double[] alpha, Vector[] polSite, double[] q, Vector[] qSite) {
            this.alpha = alpha;
            this.polSite = polSite;
            this.q = q;
            this.qSite = qSite;
        }
    }
    
    public static void main(String[] args) {
        Space space = Space3D.getInstance();
        Simulation sim = new Simulation(space);
        SpeciesSpheresRotating species = new SpeciesSpheresRotating(space, new ElementSimple("H2O", Oxygen.INSTANCE.getMass()+2*Hydrogen.INSTANCE.getMass()));
        species.setAxisSymmetric(false);
        sim.addSpecies(species);
        Box box = new Box(space);
        Box box2 = new Box(space);
        sim.addBox(box);
        box.setNMolecules(species, 3);
        box.getBoundary().setBoxSize(space.makeVector(new double[]{10000,10000,10000}));
        SpeciesWater4P water4P = new SpeciesWater4P(space);
        water4P.setConformation(new ConformationWaterGCPM(space));
        sim.addSpecies(water4P);
        sim.addBox(box2);
        box2.setNMolecules(water4P, 3);
        box2.getBoundary().setBoxSize(space.makeVector(new double[]{10000,10000,10000}));
        IAtomList triplet = box.getLeafList();
        IAtomOriented atom1 = (IAtomOriented)triplet.getAtom(0);
        IAtomOriented atom2 = (IAtomOriented)triplet.getAtom(1);
        IAtomOriented atom3 = (IAtomOriented)triplet.getAtom(2);
        AtomPair pair12 = new AtomPair(atom1, atom2);
        AtomPair pair13 = new AtomPair(atom1, atom3);
        AtomPair pair23 = new AtomPair(atom2, atom3);

        IMoleculeList moleculeTriplet = box2.getMoleculeList();
        IMolecule molecule1 = moleculeTriplet.getMolecule(0);
        IMolecule molecule2 = moleculeTriplet.getMolecule(1);
        IMolecule molecule3 = moleculeTriplet.getMolecule(2);
        MoleculePair mPair12 = new MoleculePair(molecule1, molecule2);
        MoleculePair mPair13 = new MoleculePair(molecule1, molecule3);
        MoleculePair mPair23 = new MoleculePair(molecule2, molecule3);

        P2WaterSzalewicz p2sz = new P2WaterSzalewicz(space, 2);
        p2sz.setComponent(Component.INDUCTION);
        p2sz.setBox(box);
        P2WaterSzalewicz p3sz = new P2WaterSzalewicz(space, 3);
        p3sz.setComponent(Component.INDUCTION);
        p3sz.setBox(box);

        AtomTypeAgentManager paramsManager = new AtomTypeAgentManager(null);
        P3Induction p3i = new P3Induction(space, paramsManager);
        p3i.setBox(box);
        double alphaH2O = 1.444;

        Vector polH2O = space.makeVector();
        double[] qH2O = P2WaterSzalewicz.getQ();
        Vector[] qSiteH2O = P2WaterSzalewicz.getSites(space);
        polH2O.E(qSiteH2O[0]);
        P3Induction.MyAgent agentH2O = new P3Induction.MyAgent(new double[]{alphaH2O}, new Vector[]{polH2O}, qH2O, qSiteH2O);

        paramsManager.setAgent(species.getLeafType(), agentH2O);

        PNWaterGCPM pGCPM = new PNWaterGCPM(space);
        pGCPM.setBox(box2);
        
        double r = 5;
        Vector dr1 = space.makeVector();
        Vector dr2 = space.makeVector();
        for (int i=0; i<10; i++) {
            r *= 2;
            atom2.getPosition().setX(0, r);
            atom3.getPosition().setX(0, 2*r);
            
            if (i==0) dr1.setX(0, r);
            else dr1.setX(0, r/2);
            if (i==0) dr2.setX(0, 2*r);
            else dr2.setX(0, r);
            for (int j=0; j<4; j++) {
                molecule2.getChildList().getAtom(j).getPosition().PE(dr1);
                molecule3.getChildList().getAtom(j).getPosition().PE(dr2);
            }
            
            double u2sz = p2sz.energy(pair12) + p2sz.energy(pair13) + p2sz.energy(pair23);
            double u3sz = p3sz.energy(triplet);
            double u2gcpm = pGCPM.energy(mPair12) + pGCPM.energy(mPair13) + pGCPM.energy(mPair23);
            double u3gcpm = pGCPM.energy(moleculeTriplet);
            System.out.println(r+" "+(u3sz-u2sz)+" "+p3i.energy(triplet)+" "+(u3gcpm-u2gcpm));
        }
        
        RandomMersenneTwister random = new RandomMersenneTwister(RandomNumberGeneratorUnix.getRandSeedArray());
        for (int j=0; j<10; j++) {
            System.out.println();
            r = 5;
            for (int i=0; i<10; i++) {
                r *= 2;
                atom2.getPosition().setX(0, r);
                atom3.getPosition().setX(0, r);
                atom3.getPosition().setX(1, r);
                double u2sz = p2sz.energy(pair12) + p2sz.energy(pair13) + p2sz.energy(pair23);
                double u3sz = p3sz.energy(triplet);
                double u2gcpm = pGCPM.energy(mPair12) + pGCPM.energy(mPair13) + pGCPM.energy(mPair23);
                double u3gcpm = pGCPM.energy(moleculeTriplet);
                System.out.println(r+" "+(u3sz-u2sz)+" "+p3i.energy(triplet)+" "+(u3gcpm-u2gcpm));
            }
            atom2.getOrientation().randomRotation(random, 1);
            atom3.getOrientation().randomRotation(random, 1);
        }
        
    }
}
