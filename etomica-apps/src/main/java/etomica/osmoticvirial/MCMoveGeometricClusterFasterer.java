package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.IntegratorMCFasterer;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.Potential2Soft;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.HashSet;

/**
 * Implementation of geometric cluster move of Liu and Luijten.
 * This is used to accelerate sampling of systems having molecules of very different sizes.
 * Jiwen Liu and Erik Luijten, Rejection-Free Geometric Cluster Algorithm for Complex Fluids,
 * Phys. Rev. Lett. 92, 035504, (2004).
 * DOI: 10.1103/PhysRevLett.92.035504.
 */
public class MCMoveGeometricClusterFasterer extends MCMoveBox {

    protected RandomPositionSource positionSource;
    protected AtomSource atomSource;
    protected final HashSet<IAtom> clusterAtoms;
    protected final AtomArrayList jNeighbors;
    protected final Vector oldPosition;
    protected Vector pivot;
    protected final PotentialCompute potentialCompute;
    protected final NeighborManager neighborManager;
    protected final NeighborIterator neighborIterator;
    protected final ArrayList<Vector> jNbrVecs, pairNbrVecs;
    protected final AtomArrayList atomPairs;
    protected final IntegratorMCFasterer integratorMC;
    protected final IRandom random;
    protected Potential2Soft[][] potentials;
    protected final AtomIteratorArrayListSimple atomIterator;
    protected final ISpecies solute;

    /**
     * @param space simulation space
     * @param random random number generator
     * @param integratorMC MC integrator that uses this move, needed to get the temperature
     * @param species specifies the molecules that are selected for the initial trial move; may be null, in which case
     *                any molecule in the box could be used for initial trial
     */
    public MCMoveGeometricClusterFasterer(PotentialCompute potentialCompute, NeighborManager neighborManager, Space space, IRandom random,
                                          IntegratorMCFasterer integratorMC, ISpecies species, Potential2Soft[][] pairPotentials) {
        super(null);
        clusterAtoms = new HashSet<>();
        jNeighbors = new AtomArrayList();
        this.potentialCompute = potentialCompute;
        this.neighborManager = neighborManager;
        neighborIterator = neighborManager.makeNeighborIterator();
        atomPairs = new AtomArrayList();
        jNbrVecs = new ArrayList<>();
        pairNbrVecs = new ArrayList<>();
        oldPosition = space.makeVector();
        positionSource = new RandomPositionSourceRectangular(space, random);
        if (species == null) {
            atomSource = new AtomSourceRandomLeaf();
            ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
        } else atomSource = new AtomSourceRandomSpecies(random, species);
        this.integratorMC = integratorMC;
        this.random = random;
        potentials = pairPotentials;
        atomIterator = new AtomIteratorArrayListSimple();
        this.solute = species;
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
        pivot = positionSource.randomPosition();
        IAtom atomI = atomSource.getAtom();
        if (atomI == null) return false;
        clusterAtoms.clear();
        clusterAtoms.add(atomI);
        atomPairs.clear();
        pairNbrVecs.clear();
        outer:
        while (true) {
            jNeighbors.clear();
            jNbrVecs.clear();
            gatherNeighbors(atomI);
            moveAtom(atomI.getPosition());
            neighborManager.updateAtom(atomI);
            gatherNeighbors(atomI);
            for (int x = 0; x<jNeighbors.size(); x++) {
                IAtom atomJ = jNeighbors.get(x);
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                pairNbrVecs.add(jNbrVecs.get(x));
            }
            while (atomPairs.size() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.size() - 1);
                atomI = atomPairs.remove(atomPairs.size() - 1);
                Vector drij = pairNbrVecs.remove(pairNbrVecs.size() - 1);
                if (clusterAtoms.contains(atomJ)) continue;
                double E = computeEnergy(drij, atomI, atomJ);
                moveAtom(atomJ.getPosition());
                E -= computeEnergy(drij, atomI, atomJ);
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
    private double computeEnergy(Vector dr12, IAtom atomI, IAtom atomJ) {
        int iIndex = atomI.getType().getIndex();
        int jIndex = atomJ.getType().getIndex();
        Potential2Soft p = potentials[iIndex][jIndex];
        if (p == null) {
            return 0;
        }
        return p.u(dr12, atomI, atomJ);
    }

    private void gatherNeighbors(IAtom atomI) {
        neighborIterator.iterAllNeighbors(atomI.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
            @Override
            public void accept(IAtom jAtom, Vector rij) {
                if (clusterAtoms.contains(jAtom)) return;
                jNeighbors.add(jAtom);
                Vector foo = box.getSpace().makeVector();
                foo.E(rij);
                jNbrVecs.add(foo);
            }
        });
    }

    public void setBox(Box box) {
        super.setBox(box);
        positionSource.setBox(box);
        atomSource.setBox(box);
    }

    private void moveAtom(Vector position) {
        position.TE(-1);
        position.PEa1Tv1(2, pivot);
        position.PE(box.getBoundary().centralImage(position));
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
