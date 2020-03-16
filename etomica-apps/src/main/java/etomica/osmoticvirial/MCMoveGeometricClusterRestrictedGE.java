package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.mcmove.MCMove;
import etomica.molecule.IMoleculeList;
import etomica.nbr.cell.Api1ACell;
import etomica.nbr.cell.PotentialMasterCell;
import etomica.nbr.cell.PotentialMasterCellMixed;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.Species;
import etomica.util.random.IRandom;

import java.util.*;

/**
 * Implementation of geometric cluster move for restricted gibbs ensemble by Liu, Luijten and Wilding.
 * This employs geometric cluster move (MCMoveGeometricCluster) for molecule exchange operation in a symmetrical restricted
 * Gibbs ensemble (volume of the two boxes is same and do not change over simulation).
 * Jiwen Liu, Nigel B. Wilding, and Erik Luijten, Simulation of Phase Transitions in Highly Asymmetric Fluid Mixtures,
 * Phys. Rev. Lett. 97, 115705, (2006).
 */
public class MCMoveGeometricClusterRestrictedGE extends MCMove {

    protected RandomPositionSource positionSource;
    protected AtomSource atomSource;
    protected final HashSet<IAtom> clusterAtoms1, clusterAtoms2;
    protected final List<IAtom> jNeighbors, clusterAtomsList1, clusterAtomsList2;
    protected final Vector oldPosition;
    protected Vector pivot;
    protected final Api1ACell neighbors1, neighbors2;
    protected final AtomArrayList atomPairs;
    protected final IRandom random;
    protected IPotentialAtomic[][] potentials;
    protected final AtomPair atomPair;
    protected final AtomIteratorArrayListSimple atomIterator;
    protected final Box box1, box2;
    protected double temperature;
    protected MeterPotentialEnergy energyMeter;
    protected IAtom atom;
    protected final HashMap<IAtom, Box> originalBox;
    protected final Species solute;
    protected final boolean mixedPM;

    /**
     * @param potentialMaster the PotentialMaster instance used by the simulation; this must have cell based neighbor list
     *                        which is used to identify interacting molecules
     * @param space simulation space
     * @param random random number generator
     * @param box1 specifies simulation box1
     * @param box2 specifies simulation box2
     * @param seed specifies the molecules that are selected for the initial trial move; may be null, in which case
     *                any molecule in the box could be used for initial trial
     */
    public MCMoveGeometricClusterRestrictedGE(PotentialMaster potentialMaster, Space space, IRandom random, Box box1, Box box2, Species seed) {

        super(potentialMaster);
        clusterAtoms1 = new HashSet<>();
        clusterAtoms2 = new HashSet<>();
        clusterAtomsList1 = new ArrayList<>();
        clusterAtomsList2 = new ArrayList<>();
        jNeighbors = new ArrayList<>();
        if (potentialMaster instanceof PotentialMasterCell) {
            mixedPM = potentialMaster instanceof PotentialMasterCellMixed;
            neighbors1 = new Api1ACell(((PotentialMasterCell) potentialMaster).getRange(), box1, ((PotentialMasterCell) potentialMaster).getNbrCellManager(box1));
            neighbors2 = new Api1ACell(((PotentialMasterCell) potentialMaster).getRange(), box2, ((PotentialMasterCell) potentialMaster).getNbrCellManager(box2));
        } else {
            mixedPM = false;
            neighbors1 = neighbors2 = null;
        }
        atomPairs = new AtomArrayList();
        oldPosition = space.makeVector();
        positionSource = new RandomPositionSourceRectangular(space, random);
        if(seed == null) {
            atomSource = new AtomSourceRandomLeaf();
            ((AtomSourceRandomLeaf) atomSource).setRandomNumberGenerator(random);
            System.out.println("random seeded");
        }
        else atomSource = new AtomSourceRandomSpecies(random, seed);
        this.box1 = box1;
        this.box2 = box2;
        this.random = random;
        potentials = new IPotentialAtomic[0][0];
        atomPair = new AtomPair();
        atomIterator = new AtomIteratorArrayListSimple();
        originalBox = new HashMap<>();
        energyMeter = new MeterPotentialEnergy(potentialMaster);
        temperature = 1;
        this.solute = seed;
        this.numSwaps= new int[101][2];
    }

    @Override
    public boolean doTrial() {
        int box1solute, box2solute;
        if(solute!= null) {
            box1solute = box1.getNMolecules(solute);
            box2solute = box2.getNMolecules(solute);
        }
        else{
            box1solute = box1.getMoleculeList().size();
            box2solute = box2.getMoleculeList().size();
        }
        double rand = random.nextDouble()*(box1solute+box2solute);
        Box boxI = rand < box1solute ? box1 : box2;
        atomSource.setBox(boxI);
        positionSource.setBox(boxI);
        pivot = positionSource.randomPosition();
        IAtom atomI = atomSource.getAtom();
        if (atomI == null) return false;
        clusterAtoms1.clear();
        clusterAtoms2.clear();
        clusterAtomsList1.clear();
        clusterAtomsList2.clear();
        if(boxI==box1){
            clusterAtoms2.add(atomI);
            clusterAtomsList2.add(atomI);
        }
        else{
            clusterAtoms1.add(atomI);
            clusterAtomsList1.add(atomI);
        }
        atomPairs.clear();
        outer: while(true){
            jNeighbors.clear();
            gatherNeighbors(atomI, boxI);
            originalBox.put(atomI, boxI);
            for (IAtom atomJ:jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                originalBox.put(atomJ, boxI);
            }
            jNeighbors.clear();
            moveAtom(atomI,boxI);
            Box otherBoxI = boxI==box1?box2:box1;
            gatherNeighbors(atomI, otherBoxI);
            for (IAtom atomJ:jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                originalBox.put(atomJ, otherBoxI);
            }
            while (atomPairs.size() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.size() - 1);
                atomI = atomPairs.remove(atomPairs.size() - 1);
                boxI = originalBox.get(atomI);
                Box boxJ = originalBox.get(atomJ);
                if(clusterAtoms2.contains(atomJ) || clusterAtoms1.contains(atomJ)) continue;
                double E;
                if(boxJ != boxI){
                    E = computeEnergy(atomI, atomJ);
                }
                else{
                    moveAtom(atomJ, boxJ);
                    E = -1*computeEnergy(atomI, atomJ);
                    Box otherBoxJ = boxJ==box1?box2:box1;
                    moveAtom(atomJ, otherBoxJ);
                }
                double p = 1-Math.exp(-E/temperature);
                if(p<0 || p<random.nextDouble()) continue;
                atomI = atomJ;
                boxI = boxJ;
                HashSet<IAtom> clusterAtoms = boxI==box1?clusterAtoms2:clusterAtoms1;
                clusterAtoms.add(atomI);
                List<IAtom> clusterAtomsList = boxI==box1?clusterAtomsList2:clusterAtomsList1;
                clusterAtomsList.add(atomI);
                continue outer;
            }
            break;
        }
        int newbox1solute = box1.getNMolecules(solute);
        int frac = (100 * (clusterAtomsList1.size() + clusterAtomsList2.size())) / (box1.getMoleculeList().size() + box2.getMoleculeList().size());
        if(box1solute!= newbox1solute)
            ++numSwaps[frac][1];
        ++numSwaps[frac][0];
        return true;
    }

    public int[][] numSwaps;
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
    private double computeEnergy(IAtom atomI, IAtom atomJ){
        int iIndex = atomI.getType().getIndex();
        int jIndex = atomJ.getType().getIndex();
        if(potentials.length <= iIndex){
            int oldSize = potentials.length;
            potentials = Arrays.copyOf(potentials, iIndex + 1);
            for(int i=oldSize; i<potentials.length; i++){
                potentials[i] = new IPotentialAtomic[0];
            }
        }
        if (potentials[iIndex].length <= jIndex) potentials[iIndex] = Arrays.copyOf(potentials[iIndex], jIndex + 1);
        IPotentialAtomic p = potentials[iIndex][jIndex];
        atomPair.atom0 = atomI;
        atomPair.atom1 = atomJ;
        if(p == null) {
            IPotentialAtomic ijPotential = ((PotentialMasterCell) potential).getRangedPotentials()[atomI.getType().getIndex()][atomJ.getType().getIndex()];
            potentials[iIndex][jIndex] = ijPotential;
            p = potentials[iIndex][jIndex];
            if (p == null) {
                throw new RuntimeException("could not find potential for atom types " + atomI.getType() + " " + atomJ.getType());
            }
        }
        return p.energy(atomPair);
    }

    private void gatherNeighbors(IAtom atomI, Box boxI) {
        if (neighbors1 != null && (!mixedPM || atomI.getType().getSpecies() != solute)) {
            Api1ACell neighbors = boxI == box1 ? neighbors1 : neighbors2;
            neighbors.setTarget(atomI);
            neighbors.reset();
            for (IAtomList pair = neighbors.next(); pair != null; pair = neighbors.next()) {
                IAtom atomJ = pair.get(0);
                if (atomJ == atomI) atomJ = pair.get(1);
                if (mixedPM && atomJ.getType().getSpecies() == solute) continue;
                if (clusterAtoms1.contains(atomJ) || clusterAtoms2.contains(atomJ)) continue;
                jNeighbors.add(atomJ);
            }
        }
        if (neighbors1 == null || (mixedPM && atomI.getType().getSpecies() == solute)) {
            // if mixed + solute, then loop over everything
            IAtomList allAtoms = boxI.getLeafList();
            int i = atomI.getLeafIndex();
            for (int j = 0; j < allAtoms.size(); j++) {
                if (j == i) continue;
                IAtom atomJ = allAtoms.get(j);
                if (clusterAtoms1.contains(atomJ) || clusterAtoms2.contains(atomJ)) continue;
                jNeighbors.add(atomJ);
            }
        }
        if (neighbors1 != null && mixedPM && atomI.getType().getSpecies() != solute) {
            // if mixed + solvent, loop over only solute
            IMoleculeList soluteMolecules = boxI.getMoleculeList(solute);
            for (int i = 0; i < soluteMolecules.size(); i++) {
                IAtom atomJ = soluteMolecules.get(i).getChildList().get(0);
                if (clusterAtoms1.contains(atomJ) || clusterAtoms2.contains(atomJ)) continue;
                jNeighbors.add(atomJ);
            }
        }
    }

    public void setTemperature(double temperature){
        this.temperature = temperature;
    }

    private void moveAtom(IAtom atom, Box box){
        Vector position = atom.getPosition();
        position.TE(-1);
        position.PEa1Tv1(2,pivot);
        box.removeMolecule(atom.getParentGroup());
        Box otherBox = box==box1?box2:box1;
        otherBox.addMolecule(atom.getParentGroup());
    }

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

    @Override
    public AtomIterator affectedAtoms(Box box) {
        AtomArrayList atomList = (AtomArrayList) atomIterator.getList();
        atomList.clear();
        List<IAtom> clusterAtoms = box==box1?clusterAtomsList1:clusterAtomsList2;
        for(IAtom atom: clusterAtoms){
            atomList.add(atom);
        }
        return atomIterator;
    }

    /**
     * Returns 1 because move is always accepted.
     */
    @Override
    public double energyChange(Box box) {
        return 0;
    }
}
