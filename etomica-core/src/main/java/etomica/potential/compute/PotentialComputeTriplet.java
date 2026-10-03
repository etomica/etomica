/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.exception.MethodNotImplementedException;
import etomica.integrator.IntegratorListener;
import etomica.potential.IPotential2;
import etomica.potential.IPotential3;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.SpeciesManager;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

public class PotentialComputeTriplet implements PotentialCompute {

    protected final boolean isPureAtoms;
    private final NeighborIterator neighborIterator;
    protected final IPotential3[][][] tripletPotentials;
    protected final IPotential2[][] fakePairPotentials;
    protected final Box box;
    private final NeighborManager neighborManager;
    protected Vector zero;
    protected double virialTot = Double.NaN, energyTot = Double.NaN;
    protected Vector[] forces;
    protected final Space space;
    protected final boolean centralAtom;

    /**
     * Constructs a PotentialCompute capable of iterating over triplets of atoms for 3-body potentials.
     * By default, each triplet is generated once; this works with a potential like Axilrod-Teller where each atom
     * is equivalent to the others.  With centralAtom=true, each triplet may be generated 3 times -- one
     * with each atom as a "central" atom.  This mode works with a potential like Stillinger-Weber where each
     * contribution involves a single angle.
     */
    public PotentialComputeTriplet(SpeciesManager sm, Box box, NeighborManager neighborManager, boolean centralAtom) {
        this(sm, box, neighborManager, new IPotential3[sm.getAtomTypeCount()][sm.getAtomTypeCount()][sm.getAtomTypeCount()], centralAtom);
    }

    public PotentialComputeTriplet(SpeciesManager sm, Box box, NeighborManager neighborManager, IPotential3[][][] tripletPotentials, boolean centralAtom) {
        this.centralAtom = centralAtom;
        isPureAtoms = sm.isPureAtoms();
        space = box.getSpace();
        int numAtomTypes = sm.getAtomTypeCount();
        this.neighborManager = neighborManager;
        this.tripletPotentials = tripletPotentials;
        fakePairPotentials = new IPotential2[numAtomTypes][numAtomTypes];
        for (int i=0; i<numAtomTypes; i++) {
            for (int j=i; j<numAtomTypes; j++) {
                int finalI = i, finalJ = j;
                fakePairPotentials[i][j] = fakePairPotentials[j][i] = new IPotential2() {
                    @Override
                    public double getRange() {
                        for (int k=0; k<numAtomTypes; k++) {
                            IPotential3 p3 = tripletPotentials[finalI][finalJ][k];
                            if (p3 != null) return p3.getRange();
                        }
                        return 0;
                    }
                };
            }
        }
        updateNeighborRange();


        this.neighborManager.setPairPotentials(fakePairPotentials);
        this.neighborIterator = neighborManager.makeNeighborIterator();
        this.box = box;

        zero = box.getSpace().makeVector();
        forces = new Vector[0];

    }

    @Override
    public boolean needForcesForVirial() {
        return !isPureAtoms;
    }

    public IPotential3[][][] getTripletPotentials() {
        return tripletPotentials;
    }

    @Override
    public void init() {
        this.neighborManager.init();
    }

    @Override
    public Vector[] getForces() {
        return forces;
    }

    @Override
    public double getLastVirial() {
        return virialTot;
    }

    @Override
    public double getLastEnergy() {
        return energyTot;
    }

    public void setTripletPotential(AtomType atomType1, AtomType atomType2, AtomType atomType3, IPotential3 p123) {
        tripletPotentials[atomType1.getIndex()][atomType2.getIndex()][atomType3.getIndex()] = p123;
        tripletPotentials[atomType1.getIndex()][atomType3.getIndex()][atomType2.getIndex()] = p123;
        if (!centralAtom) {
            tripletPotentials[atomType2.getIndex()][atomType1.getIndex()][atomType3.getIndex()] = p123;
            tripletPotentials[atomType2.getIndex()][atomType3.getIndex()][atomType1.getIndex()] = p123;
            tripletPotentials[atomType3.getIndex()][atomType1.getIndex()][atomType2.getIndex()] = p123;
            tripletPotentials[atomType3.getIndex()][atomType2.getIndex()][atomType1.getIndex()] = p123;
        }
        updateNeighborRange();
    }

    protected void updateNeighborRange() {
        double maxRange = Arrays.stream(fakePairPotentials).flatMap(Arrays::stream)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential2::getRange)
                .max().orElse(0);

        this.neighborManager.setPotentialRange(maxRange);
    }

    @Override
    public void updateAtom(IAtom atom) {
        this.neighborManager.updateAtom(atom);
    }

    protected final void zeroArrays(boolean doForces) {
        virialTot = 0;

        int numAtoms = box.getLeafList().size();
        if (doForces && numAtoms > forces.length) {
            int oldLength = forces.length;
            forces = Arrays.copyOf(forces, numAtoms);
            for (int i = oldLength; i < numAtoms; i++) forces[i] = box.getSpace().makeVector();
        }
        for (int i = 0; i < numAtoms; i++) {
            if (doForces) forces[i].E(0);
        }
    }

    @Override
    public double computeAll(boolean doForces, PotentialCallback pc) {
        zeroArrays(doForces);

        IAtomList atoms = box.getLeafList();
        double[] uTot = {0};
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            IPotential3[][] ip = tripletPotentials[iType];
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, (jAtom, rij, n) -> {
                int j = jAtom.getLeafIndex();
                int jType = jAtom.getType().getIndex();
                IPotential3[] pij = ip[jType];

                NeighborIterator.NeighborConsumer consumer = new NeighborIterator.NeighborConsumer() {
                    @Override
                    public void accept(IAtom kAtom, Vector rik, int m) {
                        int k = kAtom.getLeafIndex();
                        if (k<=j) return;
                        int kType = kAtom.getType().getIndex();
                        IPotential3 pijk = pij[kType];
                        if (pijk == null) return;
                        Vector rjk = space.makeVector();
                        rjk.Ev1Mv2(rij, rik);

                        double uijk;
                        double[] virial = new double[1];
                        if (doForces) {
                            uijk = pijk.udu(rij, rik, rjk, iAtom, jAtom, kAtom, virial, forces[finalI], forces[j], forces[k]);
                        }
                        else {
                            uijk = pijk.u(rij, rik, rjk, iAtom, jAtom, kAtom, virial);
                        }
                        virialTot += virial[0];
                        uTot[0] += uijk;
                    }
                };
                if (centralAtom) {
                    neighborIterator.iterAllNeighbors(finalI, consumer);
                }
                else {
                    neighborIterator.iterUpNeighbors(finalI, consumer);
                }
            });

            if (centralAtom) {
                neighborIterator.iterDownNeighbors(i, (jAtom, rij, n) -> {
                    int j = jAtom.getLeafIndex();
                    int jType = jAtom.getType().getIndex();
                    IPotential3[] pij = ip[jType];

                    NeighborIterator.NeighborConsumer consumer = new NeighborIterator.NeighborConsumer() {
                        @Override
                        public void accept(IAtom kAtom, Vector rik, int m) {
                            int k = kAtom.getLeafIndex();
                            if (k<=j) return;
                            int kType = kAtom.getType().getIndex();
                            IPotential3 pijk = pij[kType];
                            if (pijk == null) return;
                            Vector rjk = space.makeVector();
                            rjk.Ev1Mv2(rik, rij);

                            double uijk;
                            double[] virial = new double[1];
                            if (doForces) {
                                uijk = pijk.udu(rij, rik, rjk, iAtom, jAtom, kAtom, virial, forces[finalI], forces[j], forces[k]);
                            }
                            else {
                                uijk = pijk.u(rij, rik, rjk, iAtom, jAtom, kAtom, virial);
                            }
                            virialTot += virial[0];
                            uTot[0] += uijk;
                        }
                    };

                    neighborIterator.iterAllNeighbors(finalI, consumer);
                });
            }
        }

        energyTot = uTot[0];
//        System.out.println("3-body "+uTot[0]+" "+virialTot);
        return uTot[0];
    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        return computeOne(iAtom);
    }

    @Override
    public double computeOne(IAtom atom) {
        List<IAtom> seenNbrs = new ArrayList<>();
        List<Vector> seenRij = new ArrayList<>();
        return this.neighborIterator.iterAndSumAllNeighbors(atom, new NeighborIterator.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom atom2, Vector rij, int n) {
                IAtom jAtom = atom1 == atom ? atom2 : atom1;
                if (atom2==atom) rij.TE(-1);
                double uij = 0;
                for (int kk=0; kk<seenNbrs.size(); kk++) {
                    IAtom kAtom = seenNbrs.get(kk);
                    IPotential3 p3 = tripletPotentials[atom.getType().getIndex()][jAtom.getType().getIndex()][kAtom.getType().getIndex()];
                    if (p3 == null) continue;

                    Vector rjk = space.makeVector();
                    rjk.Ev1Mv2(seenRij.get(kk), rij);

                    uij += p3.u(rij, seenRij.get(kk), rjk, atom, jAtom, kAtom, new double[1]);
                }

                seenNbrs.add(jAtom);
                Vector srij = Vector.d(rij.getD());
                srij.E(rij);
                seenRij.add(srij);
                return uij;
            }
        });
    }

    @Override
    public double computeManyAtomsOld(IAtom... atoms) {
        return computeManyAtoms(atoms);
    }

    @Override
    public double computeManyAtoms(IAtom... atoms) {
        throw new MethodNotImplementedException();
    }

    @Override
    public void processAtomU(double fac) {}

    @Override
    public IntegratorListener makeIntegratorListener() {
        return this.neighborManager.makeIntegratorListener();
    }
}
