/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */
package etomica.potential.compute;

import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.box.BoxEventListener;
import etomica.box.BoxMoleculeEvent;
import etomica.integrator.IntegratorListener;
import etomica.molecule.IMolecule;
import etomica.nbr.cell.NeighborIteratorCellFasterer;
import etomica.potential.IPotential;
import etomica.potential.IPotentialEmbedding;
import etomica.potential.Potential2Soft;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.collections.DoubleArrayList;
import etomica.util.collections.IntArrayList;

import java.util.Arrays;
import java.util.Objects;

/**
 * PotentialCompute that handles EAM systems consisting of pair and embedding
 * potentials.  This class can handle purely pairwise systems, although less
 * efficiently than PotentialComputePair.
 */
public class PotentialComputeEAM implements PotentialCompute {

    private final NeighborIterator neighborIterator;
    protected final IPotentialEmbedding[] embeddingPotentials;
    protected final Potential2Soft[][] pairPotentials;
    protected final Potential2Soft[] rhoPotentials;
    protected final Box box;
    private final NeighborManager neighborManager;
    private final NeighborIteratorCellFasterer.SuperNbrConsumer nbrConsumer1, nbrConsumerEmbedding;
    protected double[] rhoSum, idf;
    protected double[] uAtom;
    protected DoubleArrayList rdrho, drhoSum;
    protected final DoubleArrayList duAtom;
    protected final IntArrayList uAtomsChanged, rhoAtomsChanged;
    protected Vector zero;
    protected double virialTot;
    protected Vector[] forces;
    protected final Space space;

    protected final int[] atomCountByType;
    protected boolean duAtomMulti = false;

    public boolean doAllTruncationCorrection = true;
    public boolean doOneTruncationCorrection = false;

    public PotentialComputeEAM(Simulation sim, Box box, NeighborManager neighborManager) {
        space = box.getSpace();
        ISpecies species = sim.getSpecies(sim.getSpeciesCount() - 1);
        int lastTypeIndex = species.getAtomType(species.getUniqueAtomTypeCount() - 1).getIndex();
        pairPotentials = new Potential2Soft[lastTypeIndex + 1][lastTypeIndex + 1];
        rhoPotentials = new Potential2Soft[lastTypeIndex + 1];
        embeddingPotentials = new IPotentialEmbedding[lastTypeIndex + 1];
        this.neighborManager = neighborManager;
        this.neighborManager.setPairPotentials(pairPotentials);
        this.neighborIterator = neighborManager.makeNeighborIterator();
        this.box = box;

        uAtom = new double[box.getLeafList().size()];
        uAtomsChanged = new IntArrayList(16);
        rhoAtomsChanged = new IntArrayList(16);
        duAtom = new DoubleArrayList(16);
        rdrho = new DoubleArrayList(16);
        drhoSum = new DoubleArrayList(16);
        zero = box.getSpace().makeVector();
        forces = new Vector[0];

        this.atomCountByType = new int[lastTypeIndex + 1];
        box.getEventManager().addListener(new BoxEventListener() {
            @Override
            public void boxMoleculeAdded(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]++;
                }
            }

            @Override
            public void boxMoleculeRemoved(BoxMoleculeEvent e) {
                for (AtomType atomType : e.getMolecule().getType().getAtomTypes()) {
                    atomCountByType[atomType.getIndex()]--;
                }
            }
        });

        // for new atom energy, includes embedding contributions
        this.nbrConsumer1 = new NeighborIteratorCellFasterer.SuperNbrConsumer() {
            @Override
            public double accept(IAtom iAtom, IAtom jAtom, Vector rij) {
                int i = iAtom.getLeafIndex();
                int j = jAtom.getLeafIndex();
                int iType = iAtom.getType().getIndex();
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = pairPotentials[iType][jType];
                double r2 = rij.squared();
                double uij = 0;
                if (pij != null) {
                    uij += pij.u(r2);
                    uAtomsChanged.add(jAtom.getLeafIndex());
                    if (duAtomMulti) {
                        duAtom.plusEquals(i, 0.5 * uij);
                        duAtom.plusEquals(j, 0.5 * uij);
                    } else {
                        duAtom.plusEquals(0, 0.5 * uij);
                        duAtom.add(0.5 * uij);
                    }
                }
                Potential2Soft irp = rhoPotentials[iType];
                Potential2Soft jrp = rhoPotentials[jType];
                if (jrp != null) {
                    double jrc = jrp.getRange();
                    if (r2 < jrc * jrc) {
                        double rho = jrp.u(r2);
                        drhoSum.plusEquals(i, rho);
                        // if i doesn't contribute to j's density, then bail
                        if (iType != jType) {
                            if (irp == null) return uij;
                            double irc = irp.getRange();
                            if (r2 > irc * irc) return uij;
                            rho = irp.u(r2);
                        }
                        rhoAtomsChanged.add(j);
                        // we need the energy of the configuration with the atom in the new spot
                        // minus the energy with no atom
                        uij += embeddingPotentials[jType].u(rhoSum[j] + drhoSum.getDouble(j) + rho);
                        uij -= embeddingPotentials[jType].u(rhoSum[j] + drhoSum.getDouble(j));
                        // drhoSum will hold new-old
                        // processAtom(+1) handles this and we should be done
                        // we'll be called again and recompute rho for the old config
                        // for processAtom(-1), we ignore drhoSum since it was already handled
                        drhoSum.plusEquals(j, rho);
                    }
                }
                return uij;
            }
        };

        // for old embedding energy
        this.nbrConsumerEmbedding = new NeighborIteratorCellFasterer.SuperNbrConsumer() {
            @Override
            public double accept(IAtom atom1, IAtom jAtom, Vector rij) {
                double r2 = rij.squared();
                int jType = jAtom.getType().getIndex();
                double jrc = rhoPotentials[jType].getRange();
                if (r2 > jrc * jrc) return 0;

                int j = jAtom.getLeafIndex();
                rhoAtomsChanged.add(j);
                double rho = rhoPotentials[jType].u(r2);
                drhoSum.replace(j, -rho);
                // we need to compute the old energy of the system with and without iAtom
                return embeddingPotentials[jType].u(rhoSum[j])
                        - embeddingPotentials[jType].u(rhoSum[j] - rho);
            }
        };
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
    public double getOldEnergy() {
        double uTot = 0;
        for (double iuAtom : uAtom) {
            uTot += iuAtom;
        }
        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot += uCorrection[0];
        return uTot;
    }

    public void setPairPotential(AtomType atomType1, AtomType atomType2, Potential2Soft p12) {
        pairPotentials[atomType1.getIndex()][atomType2.getIndex()] = p12;
        pairPotentials[atomType2.getIndex()][atomType1.getIndex()] = p12;

        updateMaxRange();
    }

    public void setRhoPotential(AtomType atomType, Potential2Soft pRho) {
        rhoPotentials[atomType.getIndex()] = pRho;

        updateMaxRange();
    }

    private void updateMaxRange() {

        double maxRange1 = Arrays.stream(pairPotentials).flatMap(Arrays::stream)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential::getRange)
                .max().orElse(0);

        double maxRange2 = Arrays.stream(rhoPotentials)
                .filter(Objects::nonNull)
                .mapToDouble(IPotential::getRange)
                .max().orElse(0);

        this.neighborManager.setPotentialRange(Math.max(maxRange1, maxRange2));
    }

    public void setEmbeddingPotential(AtomType atomType, IPotentialEmbedding p) {
        embeddingPotentials[atomType.getIndex()] = p;
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
        if (numAtoms > uAtom.length) {
            uAtom = new double[numAtoms];
            rhoSum = new double[numAtoms];
            idf = new double[numAtoms];
        }
        for (int i = 0; i < numAtoms; i++) {
            uAtom[i] = 0;
            rhoSum[i] = 0;
            if (doForces) forces[i].E(0);
        }
    }

    @Override
    public double computeAll(boolean doForces) {
        zeroArrays(doForces);
        rdrho.clear();

        IAtomList atoms = box.getLeafList();
        double[] uTot = {0};
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            Potential2Soft[] ip = pairPotentials[iType];
            Potential2Soft irp = rhoPotentials[iType];
            int finalI = i;
            neighborIterator.iterUpNeighbors(i, (jAtom, rij) -> {
                int j = jAtom.getLeafIndex();
                int jType = jAtom.getType().getIndex();
                Potential2Soft pij = ip[jType];
                double[] u012 = new double[3];
                double r2 = rij.squared();
                if (pij != null) {
                    pij.u012add(r2, u012);
                    double uij = u012[0];
                    if (uij != 0) {
                        uAtom[finalI] += 0.5 * uij;
                        uAtom[j] += 0.5 * uij;
                        double duij = u012[1];
                        virialTot += duij;
                        if (doForces) {
                            rij.TE(duij / r2);
                            forces[finalI].PE(rij);
                            forces[j].ME(rij);
                        }
                        uTot[0] += uij;
                    }

                    u012[0] = u012[1] = u012[2];
                }
                if (irp != null && embeddingPotentials[jType] != null) {
                    // i contributes to j density
                    double irc = irp.getRange();
                    u012[0] = 0;
                    if (r2 < irc * irc) {
                        // i contributes to j density
                        irp.u012add(r2, u012);
                        rhoSum[j] += u012[0];
                        rdrho.add(u012[1]);
                    }

                    if (iType != jType) {
                        Potential2Soft jrp = rhoPotentials[jType];
                        if (jrp != null && embeddingPotentials[iType] != null) {
                            double jrc = jrp.getRange();
                            if (r2 < jrc * jrc) {
                                // j contributes to i density
                                jrp.u012add(r2, u012);
                                rdrho.add(u012[1]);
                                rhoSum[finalI] += u012[0];
                            }
                        }
                    } else {
                        rhoSum[finalI] += u012[0];
                    }
                } else if (iType != jType) {
                    Potential2Soft jrp = rhoPotentials[jType];
                    if (jrp != null && embeddingPotentials[iType] != null) {
                        double jrc = jrp.getRange();
                        if (r2 < jrc * jrc) {
                            // j contributes to i density
                            jrp.u012add(r2, u012);
                            rdrho.add(u012[1]);
                            rhoSum[finalI] += u012[0];
                        }
                    }
                }
            });
        }

        int[] rdrhoIdx = {0};
        for (int i = 0; i < atoms.size(); i++) {
            IAtom iAtom = atoms.get(i);
            int iType = iAtom.getType().getIndex();
            if (embeddingPotentials[iType] == null) continue;
            double[] u = {0}, du = {0};
            embeddingPotentials[iType].udu(rhoSum[i], u, du);
            uTot[0] += u[0];
            if (doForces) {
                idf[i] = du[0];
            }
        }
        if (doForces) {
            for (int i = 0; i < atoms.size(); i++) {
                IAtom iAtom = atoms.get(i);
                int iType = iAtom.getType().getIndex();

                Potential2Soft irp = rhoPotentials[iType];
                double iCutoff2 = irp.getRange() * irp.getRange();
                int finalI = i;
                neighborIterator.iterUpNeighbors(i, (jAtom, rij) -> {
                    int j = jAtom.getLeafIndex();
                    int jType = jAtom.getType().getIndex();
                    if (rhoPotentials[jType] == null) return;
                    double r2 = rij.squared();

                    if (r2 < iCutoff2) {
                        double fac;
                        if (iType == jType) {
                            fac = (idf[finalI] + idf[j]) * rdrho.getDouble(rdrhoIdx[0]);
                            rdrhoIdx[0]++;
                        } else {
                            fac = idf[finalI] * rdrho.getDouble(rdrhoIdx[0]);
                            rdrhoIdx[0]++;
                            double jrc = rhoPotentials[jType].getRange();
                            if (r2 < jrc * jrc) {
                                fac += idf[j] * rdrho.getDouble(rdrhoIdx[0]);
                                rdrhoIdx[0]++;
                            }
                        }
                        virialTot += fac;
                        fac /= r2;
                        rij.TE(fac);
                        forces[finalI].PE(rij);
                        forces[j].ME(rij);
                    }
                });
            }
        }

        double[] uCorrection = new double[1];
        double[] duCorrection = new double[1];
        this.computeAllTruncationCorrection(uCorrection, duCorrection);
        uTot[0] += uCorrection[0];
        virialTot += duCorrection[0];

        return uTot[0];
    }

    public double oldEmbeddingEnergy(IAtom iAtom) {
        // just compute all the embedding energies
        int i = iAtom.getLeafIndex();
        rhoAtomsChanged.clear();
        rhoAtomsChanged.add(i);
        int iType = iAtom.getType().getIndex();
        double u = embeddingPotentials[iType].u(rhoSum[i]);
        u += this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumerEmbedding);
        return u;

    }

    @Override
    public double computeOneOld(IAtom iAtom) {
        double u = this.computeOneTruncationCorrection(iAtom.getLeafIndex());
        return u + uAtom[iAtom.getLeafIndex()] * 2 + oldEmbeddingEnergy(iAtom);
    }

    @Override
    public double computeOneOldMolecule(IMolecule molecule) {
        double u = 0;
        for (IAtom atom : molecule.getChildList()) {
            u += uAtom[atom.getLeafIndex()] * 2;
            u += this.computeOneTruncationCorrection(atom.getLeafIndex());
        }
        return u;
    }

    @Override
    public double computeOne(IAtom iAtom) {
        this.duAtomMulti = false;
        int i = iAtom.getLeafIndex();
        IAtomList atoms = box.getLeafList();
        uAtomsChanged.clear();
        uAtomsChanged.ensureCapacity(atoms.size());
        uAtomsChanged.add(i);
        duAtom.clear();
        duAtom.ensureCapacity(atoms.size());
        duAtom.add(0);
        drhoSum.clear();
        drhoSum.ensureCapacity(atoms.size());
        if (rhoAtomsChanged.size() == 0) {
            // if we called oldEnergy for this atom, then it will already be in the list
            rhoAtomsChanged.add(i);
        }
        double u = this.computeOneInternal(iAtom);
        u += this.computeOneTruncationCorrection(i);
        return u;
    }

    protected double computeOneInternal(IAtom iAtom) {
        return this.neighborIterator.iterAndSumAllNeighbors(iAtom, this.nbrConsumer1);
    }

    @Override
    public double computeOneMolecule(IMolecule molecule) {
        duAtomMulti = true;
        uAtomsChanged.clear();
        duAtom.setAll(0);
        duAtom.ensureCapacity(box.getLeafList().size());
        double u = 0;

        for (IAtom atom : molecule.getChildList()) {
            uAtomsChanged.add(atom.getLeafIndex());

            u += computeOneInternal(atom);
            u += computeOneTruncationCorrection(atom.getLeafIndex());
        }

        return u;
    }

    @Override
    public void processAtomU(double fac) {
        for (int j = 0; j < uAtomsChanged.size(); j++) {
            int jj = uAtomsChanged.getInt(j);
            if (duAtomMulti) {
                uAtom[jj] += fac * duAtom.getDouble(jj);
                duAtom.replace(jj, 0);
            } else {
                uAtom[jj] += fac * duAtom.getDouble(j);
            }
        }

        for (int j = 0; j < rhoAtomsChanged.size(); j++) {
            int jj = rhoAtomsChanged.getInt(j);
            // by the time we get to processAtom(+1), we have
            // the difference.  just use that
            if (fac == 1) {
                rhoSum[jj] += drhoSum.getDouble(jj);
            }
        }
    }

    @Override
    public IntegratorListener makeIntegratorListener() {
        return this.neighborManager.makeIntegratorListener();
    }

    public void computeAllTruncationCorrection(double[] uCorrection, double[] duCorrection) {
        if (!doAllTruncationCorrection) {
            return;
        }
        double uCor = 0;
        double duCor = 0;
        for (int i = 0; i < atomCountByType.length; i++) {
            for (int j = i; j < atomCountByType.length; j++) {
                Potential2Soft p = pairPotentials[i][j];
                int numPairs;
                if (j == i) {
                    numPairs = atomCountByType[i] * (atomCountByType[j] - 1) / 2;
                } else {
                    numPairs = atomCountByType[i] * atomCountByType[j];
                }
                double pairDensity = numPairs / box.getBoundary().volume();

                double[] u = new double[1];
                double[] du = new double[1];
                p.u01TruncationCorrection(u, du);
                uCor += pairDensity * u[0];
                duCor += pairDensity * du[0];

            }
        }
        uCorrection[0] = uCor;
        duCorrection[0] = duCor;
    }

    public double computeOneTruncationCorrection(int iAtom) {
        if (!doOneTruncationCorrection) {
            return 0;
        }
        int iType = box.getLeafList().get(iAtom).getType().getIndex();
        double uCorrection = 0;
        for (int j = 0; j < atomCountByType.length; j++) {
            Potential2Soft p = pairPotentials[iType][j];
            double pairDensity;
            if (iType == j) {
                pairDensity = (atomCountByType[j] - 1) / box.getBoundary().volume();
            } else {
                pairDensity = atomCountByType[j] / box.getBoundary().volume();
            }
            double integral = p.integral(p.getRange());
            uCorrection += pairDensity * integral;
        }
        return uCorrection;
    }
}
