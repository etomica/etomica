package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayList;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.nbr.NeighborCriterion;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.potential.IPotential;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialArray;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.Arrays;
import etomica.util.random.IRandom;

import java.util.HashSet;
import java.util.Random;

/**
 * Implementation of geometric cluster move of Liu and Luijten.
 * This is used to accelerate sampling of systems having molecules of very different sizes.
 * Jiwen Liu and Erik Luijten, Rejection-Free Geometric Cluster Algorithm for Complex Fluids,
 * Phys. Rev. Lett. 92, 035504, (2004).
 * DOI: 10.1103/PhysRevLett.92.035504.
 */
public class MCMoveGeometricCluster extends MCMoveBox {

    protected RandomPositionSource positionSource;
    protected AtomSource atomSource;
    protected final HashSet<IAtom> clusterAtoms, jNeighbors;
    protected final Vector oldPosition;
    protected Vector pivot;
    protected final Api1ACell neighbors;
    protected final AtomArrayList atomPairs;
    protected final IntegratorMC integratorMC;
    protected final IRandom random;
    protected IPotentialAtomic[][] potentials;
    protected final AtomPair atomPair;
    protected final AtomIteratorArrayListSimple atomIterator;

    /**
     * @param potentialMaster the PotentialMaster instance used by the simulation; this must have cell based neighbor list
     *                        which is used to identify interacting molecules
     * @param space simulation space
     * @param random random number generator
     * @param integratorMC MC integrator that uses this move, needed to get the temperature
     * @param species specifies the molecules that are selected for the initial trial move; may be null, in which case
     *                any molecule in the box could be used for initial trial
     */
    public MCMoveGeometricCluster(PotentialMasterCell potentialMaster, Space space, IRandom random
            , IntegratorMC integratorMC, ISpecies species) {
        super(potentialMaster);
        clusterAtoms = new HashSet<>();
        jNeighbors = new HashSet<>();
        neighbors = new Api1ACell(space.D(), potentialMaster.getRange(), potentialMaster.getCellAgentManager());
        atomPairs = new AtomArrayList();
        oldPosition = space.makeVector();
        positionSource = new RandomPositionSourceRectangular(space, random);
        if (species == null) {
            atomSource = new AtomSourceRandomLeaf();
            ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
        } else atomSource = new AtomSourceRandomSpecies(random, species);
        neighbors.setDirection(null);
        this.integratorMC = integratorMC;
        this.random = random;
        potentials = new IPotentialAtomic[0][0];
        atomPair = new AtomPair();
        atomIterator = new AtomIteratorArrayListSimple();
    }

    @Override
    public AtomIterator affectedAtoms() {
        AtomArrayList atomList = (AtomArrayList) atomIterator.getList();
        atomList.clear();
        for (IAtom atom : clusterAtoms) {
            atomList.add(atom);
        }
        return atomIterator;
    }

    @Override
    public double energyChange() {
        return 0;
    } //this isn't correct in general

    @Override
    public boolean doTrial() {
        NeighborCellManager ncm = ((PotentialMasterCell) potential).getNbrCellManager(box);
        pivot = positionSource.randomPosition();
        IAtom atomI = atomSource.getAtom();
        if (atomI == null) return false;
        clusterAtoms.clear();
        clusterAtoms.add(atomI);
        atomPairs.clear();
        outer:
        while (true) {
            jNeighbors.clear();
            gatherNeighbors(atomI);
            moveAtom(atomI.getPosition());
            ncm.updateCell(atomI);
            gatherNeighbors(atomI);
            for (IAtom atomJ : jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
            }
            while (atomPairs.getAtomCount() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.getAtomCount() - 1);
                atomI = atomPairs.remove(atomPairs.getAtomCount() - 1);
                if (clusterAtoms.contains(atomJ)) continue;
                double E = computeEnergy(atomI, atomJ);
                moveAtom(atomJ.getPosition());
                E -= computeEnergy(atomI, atomJ);
                double p = 1 - Math.exp(-E / integratorMC.getTemperature());
                moveAtom(atomJ.getPosition());
                if (p < 0 || p < random.nextDouble()) continue;
                atomI = atomJ;
                clusterAtoms.add(atomI);
                continue outer;
            }
            break;
        }
        return true;
    }

   //Computes the energy between atomI and atomJ, determines the potential as needed and saves for later use.
    private double computeEnergy(IAtom atomI, IAtom atomJ) {
        int iIndex = atomI.getType().getIndex();
        int jIndex = atomJ.getType().getIndex();
        if (potentials.length <= iIndex) {
            int oldSize = potentials.length;
            potentials = (IPotentialAtomic[][]) Arrays.resizeArray(potentials, iIndex + 1);
            for (int i = oldSize; i < potentials.length; i++) {
                potentials[i] = new IPotentialAtomic[0];
            }
        }
        if (potentials[iIndex].length <= jIndex)
            potentials[iIndex] = (IPotentialAtomic[]) Arrays.resizeArray(potentials[iIndex], jIndex + 1);
        IPotentialAtomic p = potentials[iIndex][jIndex];
        atomPair.atom0 = atomI;
        atomPair.atom1 = atomJ;
        if (p == null) {
            PotentialArray iPotentials = ((PotentialMasterCell) potential).getRangedPotentials(atomI.getType());
            NeighborCriterion[] neighborCriteria = iPotentials.getCriteria();
            IPotential[] iPotential = iPotentials.getPotentials();
            for (int i = 0; i < iPotential.length; i++) {
                if (neighborCriteria[i].accept(atomPair)) {
                    potentials[iIndex][jIndex] = (IPotentialAtomic) iPotential[i];
                    p = potentials[iIndex][jIndex];
                }
            }
        }
        return p.energy(atomPair);
    }

    private void gatherNeighbors(IAtom atomI) {
        neighbors.setTarget(atomI);
        neighbors.reset();
        for (IAtomList pair = neighbors.next(); pair != null; pair = neighbors.next()) {
            IAtom atomJ = pair.getAtom(0);
            if (atomJ == atomI) atomJ = pair.getAtom(1);
            if (clusterAtoms.contains(atomJ)) continue;
            jNeighbors.add(atomJ);
        }
    }

    public void setBox(Box box) {
        super.setBox(box);
        neighbors.setBox(box);
        positionSource.setBox(box);
        atomSource.setBox(box);
    }

    private void moveAtom(Vector position) {
        position.TE(-1);
        position.PEa1Tv1(2, pivot);
    }

    /**
     * Returns 1 because move is always accepted.
     */
    @Override
    public double getChi(double temperature) {
        return 1;
    }

    @Override
    public void acceptNotify() {

    }

    @Override
    public void rejectNotify() {

    }
}
