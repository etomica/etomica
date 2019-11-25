package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.NeighborCellManager;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.cell.PotentialMasterCellMixed;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

import java.util.Arrays;
import java.util.HashSet;

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
    protected Api1ACell neighbors;
    protected final AtomArrayList atomPairs;
    protected final IntegratorMC integratorMC;
    protected final IRandom random;
    protected IPotentialAtomic[][] potentials;
    protected final AtomPair atomPair;
    protected final AtomIteratorArrayListSimple atomIterator;
    protected final ISpecies solute;
    protected final boolean mixedPM;

    /**
     * @param potentialMaster the PotentialMaster instance used by the simulation; this must have cell based neighbor list
     *                        which is used to identify interacting molecules
     * @param space simulation space
     * @param random random number generator
     * @param integratorMC MC integrator that uses this move, needed to get the temperature
     * @param species specifies the molecules that are selected for the initial trial move; may be null, in which case
     *                any molecule in the box could be used for initial trial
     */
    public MCMoveGeometricCluster(PotentialMaster potentialMaster, Space space, IRandom random
            , IntegratorMC integratorMC, ISpecies species) {
        super(potentialMaster);
        clusterAtoms = new HashSet<>();
        jNeighbors = new HashSet<>();
        if (potentialMaster instanceof PotentialMasterCell) {
            mixedPM = potentialMaster instanceof PotentialMasterCellMixed;
        } else {
            mixedPM = false;
        }
        atomPairs = new AtomArrayList();
        oldPosition = space.makeVector();
        positionSource = new RandomPositionSourceRectangular(space, random);
        if (species == null) {
            atomSource = new AtomSourceRandomLeaf();
            ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
        } else atomSource = new AtomSourceRandomSpecies(random, species);
        this.integratorMC = integratorMC;
        this.random = random;
        potentials = new IPotentialAtomic[0][0];
        atomPair = new AtomPair();
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
            while (atomPairs.size() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.size() - 1);
                atomI = atomPairs.remove(atomPairs.size() - 1);
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

    /**
     * Informs the move that p2 is the potential between type1 and type2.
     * type1 and type2 may the same.
     *
     * @param type1 the first atom type
     * @param type2 the second atom type
     * @param p2    the potential between type1 and type2
     */
    public void setPotential(AtomType type1, AtomType type2, IPotentialAtomic p2) {
        int iIndex = type1.getIndex();
        int jIndex = type2.getIndex();
        int maxIndex = iIndex > jIndex ? iIndex : jIndex;
        if (potentials.length <= maxIndex) {
            int oldSize = potentials.length;
            potentials = Arrays.copyOf(potentials, maxIndex + 1);
            for (int i = 0; i < oldSize; i++) {
                potentials[i] = Arrays.copyOf(potentials[i], maxIndex + 1);
            }
            for (int i = oldSize; i < potentials.length; i++) {
                potentials[i] = new IPotentialAtomic[maxIndex + 1];
            }
        }
        potentials[iIndex][jIndex] = potentials[jIndex][iIndex] = p2;
    }

    //Computes the energy between atomI and atomJ, determines the potential as needed and saves for later use.
    private double computeEnergy(IAtom atomI, IAtom atomJ) {
        int iIndex = atomI.getType().getIndex();
        int jIndex = atomJ.getType().getIndex();
        if (potentials.length <= iIndex) {
            int oldSize = potentials.length;
            potentials = Arrays.copyOf(potentials, iIndex + 1);
            for (int i = oldSize; i < potentials.length; i++) {
                potentials[i] = new IPotentialAtomic[0];
            }
        }
        if (potentials[iIndex].length <= jIndex)
            potentials[iIndex] = Arrays.copyOf(potentials[iIndex], jIndex + 1);
        IPotentialAtomic p = potentials[iIndex][jIndex];
        if (p == null) {
            IPotentialAtomic ijPotential = ((PotentialMasterCell) potential).getRangedPotentials()[atomI.getType().getIndex()][atomJ.getType().getIndex()];
            potentials[iIndex][jIndex] = ijPotential;
            p = potentials[iIndex][jIndex];
            if (p == null) {
                throw new RuntimeException("could not find potential for atom types " + atomI.getType() + " " + atomJ.getType());
            }
        }
        atomPair.atom0 = atomI;
        atomPair.atom1 = atomJ;
        return p.energy(atomPair);
    }

    private void gatherNeighbors(IAtom atomI) {
        if (neighbors != null && (!mixedPM || atomI.getType().getSpecies() != solute)) {
            neighbors.setTarget(atomI);
            neighbors.reset();
            for (IAtomList pair = neighbors.next(); pair != null; pair = neighbors.next()) {
                IAtom atomJ = pair.get(0);
                if (atomJ == atomI) atomJ = pair.get(1);
                if (mixedPM && atomJ.getType().getSpecies() == solute) continue;
                if (clusterAtoms.contains(atomJ)) continue;
                jNeighbors.add(atomJ);
            }
        }
        if (neighbors == null || (mixedPM && atomI.getType().getSpecies() == solute)) {
            // if mixed + solute, then loop over everything
            IAtomList allAtoms = box.getLeafList();
            int i = atomI.getLeafIndex();
            for (int j = 0; j < allAtoms.size(); j++) {
                if (j == i) continue;
                IAtom atomJ = allAtoms.get(j);
                if (clusterAtoms.contains(atomJ)) continue;
                jNeighbors.add(atomJ);
            }
        }
        if (neighbors != null && mixedPM && atomI.getType().getSpecies() != solute) {
            // if mixed + solvent, loop over only solute
            IMoleculeList soluteMolecules = box.getMoleculeList(solute);
            for (int i = 0; i < soluteMolecules.size(); i++) {
                IAtom atomJ = soluteMolecules.get(i).getChildList().get(0);
                if (clusterAtoms.contains(atomJ)) continue;
                jNeighbors.add(atomJ);
            }
        }
    }

    public void setBox(Box box) {
        super.setBox(box);
        if (potential instanceof PotentialMasterCell) {
            neighbors = new Api1ACell(((PotentialMasterCell) potential).getRange(), box, ((PotentialMasterCell) potential).getNbrCellManager(box));
            neighbors.setDirection(null);
        }
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
