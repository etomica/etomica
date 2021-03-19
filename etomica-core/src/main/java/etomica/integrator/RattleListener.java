/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.integrator;

import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.ShakeListener.BondConstraints;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Boundary;
import etomica.space.Vector;
import etomica.species.SpeciesAgentManager;
import etomica.species.SpeciesManager;

import java.util.Arrays;

public class RattleListener implements IntegratorListenerMD {

    protected final IntegratorMDFasterer integrator;
    protected Vector[][] drOld;
    protected Vector[] rOld;
    protected final SpeciesAgentManager<BondConstraints> agentManager;
    protected double shakeTol;
    protected int maxIterations;
    protected final int[] nBonds;
    protected boolean[][] moved;
    protected boolean ready = false;

    public RattleListener(SpeciesManager sm, SpeciesAgentManager<BondConstraints> agentManager, IntegratorMDFasterer integrator) {
        this.integrator = integrator;
        this.agentManager = agentManager;
        nBonds = new int[sm.getSpeciesCount()];
        for (int i = 0; i < nBonds.length; i++) {
            BondConstraints bc = agentManager.getAgent(sm.getSpecies(i));
            nBonds[i] = bc == null ? 0 : bc.bondedAtoms.length;
        }
        moved = new boolean[2][0];
        drOld = new Vector[0][];
        rOld = new Vector[0];
        setShakeTolerance(1e-12);
        setMaxIterations(20);
    }

    public void setShakeTolerance(double newShakeTol) {
        shakeTol = newShakeTol;
    }

    public double getShakeTolerance() {
        return shakeTol;
    }

    public int getMaxIterations() {
        return maxIterations;
    }

    public void setMaxIterations(int newMaxIterations) {
        maxIterations = newMaxIterations;
    }

    public void integratorStepStarted(IntegratorEvent e) {
        Box box = integrator.getBox();
        IMoleculeList molecules = box.getMoleculeList();
        if (drOld.length < molecules.size()) {
            int oldLength = drOld.length;
            drOld = Arrays.copyOf(drOld, molecules.size());
            for (int i = oldLength; i < drOld.length; i++) {
                drOld[i] = new Vector[0];
            }
        }
        for (int i = 0; i < molecules.size(); i++) {
            IMolecule m = molecules.get(i);
            int nb = nBonds[m.getType().getIndex()];
            if (drOld[i].length < nb) {
                int oldLength = drOld[i].length;
                drOld[i] = Arrays.copyOf(drOld[i], nb);
                for (int j = oldLength; j < nb; j++) {
                    drOld[i][j] = box.getSpace().makeVector();
                }
            }
        }
        if (rOld.length < box.getLeafList().size()) {
            int oldLength = rOld.length;
            rOld = Arrays.copyOf(rOld, box.getLeafList().size());
            for (int i = oldLength; i < rOld.length; i++) {
                rOld[i] = box.getSpace().makeVector();
            }
        }
        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            BondConstraints bondConstraints = agentManager.getAgent(molecule.getType());
            if (bondConstraints == null) continue;
            IAtomList childList = molecule.getChildList();
            Boundary boundary = box.getBoundary();

            for (int j = 0; j < bondConstraints.bondedAtoms.length; j++) {
                IAtom atom0 = childList.get(bondConstraints.bondedAtoms[j][0]);
                IAtom atom1 = childList.get(bondConstraints.bondedAtoms[j][1]);
                drOld[i][j].Ev1Mv2(atom1.getPosition(), atom0.getPosition());
                boundary.nearestImage(drOld[i][j]);
                double bl = Math.sqrt(drOld[i][j].squared());
//                System.out.printf("%3d %1d %+6.1e %6.3f %6.3f\n", i, j, bl-bondConstraints.bondLengths[j], bl, bondConstraints.bondLengths[j]);
//                if (Math.abs(bl-bondConstraints.bondLengths[j]) > 1e-11) throw new RuntimeException("oops");
            }
        }
        for (IAtom a : box.getLeafList()) {
            rOld[a.getLeafIndex()].E(a.getPosition());
        }
        ready = true;
    }


    public void integratorForcePrecomputed(IntegratorEvent e) {
        if (!ready) return;
        Box box = integrator.getBox();
        IMoleculeList molecules = box.getMoleculeList();
        Boundary boundary = box.getBoundary();

        Vector dr = box.getSpace().makeVector();
        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            BondConstraints bondConstraints = agentManager.getAgent(molecule.getType());
            if (bondConstraints == null) continue;
            IAtomList childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            double[] bondLengths = bondConstraints.bondLengths;

            if (childList.size() > moved[0].length) {
                moved = new boolean[2][childList.size()];
            }
            for (int j = 0; j < childList.size(); j++) {
                moved[1][j] = true;
            }

            for (int iter = 0; iter < maxIterations; iter++) {
                boolean success = true;
                for (int j = 0; j < childList.size(); j++) {
                    moved[0][j] = moved[1][j];
                    moved[1][j] = false;
                }
                double maxDiff = 0;
                for (int j = 0; j < bondConstraints.bondedAtoms.length; j++) {
                    int iAtom1 = bondedAtoms[j][0];
                    int iAtom2 = bondedAtoms[j][1];
                    if (!moved[0][iAtom1] && !moved[0][iAtom2]) {
                        continue;
                    }
                    IAtomKinetic atom1 = (IAtomKinetic) childList.get(iAtom1);
                    IAtomKinetic atom2 = (IAtomKinetic) childList.get(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
                    //                    if (i==0) System.out.println(iter+" old dr "+Math.sqrt(dr.squared())+" vs "+bondLengths[j]);
                    double dr2 = dr.squared();
                    double bl2 = bondLengths[j] * bondLengths[j];
                    double diffSq = bl2 - dr2;
                    if (diffSq / bl2 > maxDiff) maxDiff = diffSq / bl2;
                    if (Math.abs(diffSq / bl2) > shakeTol) {
                        double mass1 = atom1.getType().getMass();
                        double mass2 = atom2.getType().getMass();
                        double rMass = 1.0 / mass1 + 1.0 / mass2;
                        double drDotDrOld = dr.dot(drOld[i][j]);
                        if (drDotDrOld / bl2 < 0.1) {
                            System.out.println("molecule " + i + " bond " + j);
                            System.out.println("atoms " + atom1 + " " + atom1.getType() + " " + atom2 + " " + atom2.getType());
                            System.out.println("dr " + dr);
                            System.out.println("drOld " + drOld[i][j]);
                            System.out.println("drDotDrOld " + drDotDrOld);
                            throw new RuntimeException("oops");
                        }
                        double gab = diffSq / (2.0 * rMass * drDotDrOld);
                        atom2.getPosition().PEa1Tv1(gab / mass2, drOld[i][j]);
                        atom1.getPosition().PEa1Tv1(-gab / mass1, drOld[i][j]);

                        gab /= integrator.getTimeStep();

                        atom2.getVelocity().PEa1Tv1(gab / mass2, drOld[i][j]);
                        atom1.getVelocity().PEa1Tv1(-gab / mass1, drOld[i][j]);

                        moved[1][iAtom1] = true;
                        moved[1][iAtom2] = true;
                        success = false;
                    }
                }
                if (success) {
                    break;
                }
//                System.out.println(i+" "+iter+" "+maxDiff);
                if (iter == maxIterations - 1) {
                    System.err.println("failed to converge in shake for molecule " + i + " " + maxDiff);
                }
            }
        }

        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            BondConstraints bondConstraints = agentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            bondConstraints.relaxMolecule(molecule);
        }

    }

    @Override
    public void integratorForceComputed(IntegratorEvent e) {
        if (!ready) return;
        Box box = integrator.getBox();
        IMoleculeList molecules = box.getMoleculeList();
        Vector[] forces = integrator.getPotentialCompute().getForces();
        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            BondConstraints bondConstraints = agentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }
            bondConstraints.redistributeForces(molecule, forces);
        }
    }

    @Override
    public void preThermostat(IntegratorEvent e) {
        /*
         * Rattle Part II
         */
        if (!ready) return;
        int numIterations = 0;
        Box box = integrator.getBox();
        IMoleculeList molecules = box.getMoleculeList();
        Vector dr = box.getSpace().makeVector();
        Vector dv = box.getSpace().makeVector();
        for (int i = 0; i < molecules.size(); i++) {
            IMolecule molecule = molecules.get(i);
            BondConstraints bondConstraints = agentManager.getAgent(molecule.getType());
            if (bondConstraints == null) {
                continue;
            }

            IAtomList childList = molecule.getChildList();
            int[][] bondedAtoms = bondConstraints.bondedAtoms;
            Boundary boundary = box.getBoundary();
            double[] bondLengths = bondConstraints.bondLengths;

            for (int j = 0; j < childList.size(); j++) {
                moved[1][j] = true;
            }

            for (int iter = 0; iter < maxIterations; iter++) {
                numIterations++;
                boolean success = true;
                for (int j = 0; j < childList.size(); j++) {
                    moved[0][j] = moved[1][j];
                    moved[1][j] = false;
                }
                for (int j = 0; j < bondConstraints.bondedAtoms.length; j++) {
                    int iAtom1 = bondedAtoms[j][0];
                    int iAtom2 = bondedAtoms[j][1];
                    if (!moved[0][iAtom1] && !moved[0][iAtom2]) {
                        continue;
                    }
                    IAtomKinetic atom1 = (IAtomKinetic) childList.get(iAtom1);
                    IAtomKinetic atom2 = (IAtomKinetic) childList.get(iAtom2);
                    dr.Ev1Mv2(atom2.getPosition(), atom1.getPosition());
                    boundary.nearestImage(dr);
                    dv.Ev1Mv2(atom2.getVelocity(), atom1.getVelocity());
                    double drdotdv = dr.dot(dv);
                    double mass1 = atom1.getType().getMass();
                    double mass2 = atom2.getType().getMass();
                    double bl2 = bondLengths[j] * bondLengths[j];
                    double g = -drdotdv / ((1.0 / mass1 + 1.0 / mass2) * bl2);
                    if (Math.abs(g) > shakeTol) {
                        dr.TE(g);

                        atom2.getVelocity().PEa1Tv1(1.0 / (mass2), dr);
                        atom1.getVelocity().PEa1Tv1(-1.0 / (mass1), dr);

                        moved[1][iAtom1] = true;
                        moved[1][iAtom2] = true;
                        success = false;
                    }
                }
                if (success) {
                    break;
                }
                if (iter == maxIterations - 1) {
                    System.err.println("failed to converge in rattle for molecule " + i);
                }
            }
        }

    }
}
