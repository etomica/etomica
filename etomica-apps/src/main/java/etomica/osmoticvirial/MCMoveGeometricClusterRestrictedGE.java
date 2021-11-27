package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.mcmove.MCMove;
import etomica.potential.IPotential2;
import etomica.potential.compute.NeighborIterator;
import etomica.potential.compute.NeighborManager;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.species.ISpecies;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

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
    protected final NeighborManager neighborManager1, neighborManager2;
    protected final NeighborIterator neighborIterator1, neighborIterator2;
    protected final ArrayList<Vector> jNbrVecs, pairNbrVecs;
    protected final AtomArrayList atomPairs;
    protected final IRandom random;
    protected IPotential2[][] potentials;
    protected final Box box1, box2;
    protected double temperature;
    protected IAtom atom;
    protected final HashMap<IAtom, Box> originalBox;
    protected final ISpecies solute;

    /**
     * @param space simulation space
     * @param random random number generator
     * @param box1 specifies simulation box1
     * @param box2 specifies simulation box2
     * @param seed specifies the molecules that are selected for the initial trial move; may be null, in which case
     *                any molecule in the box could be used for initial trial
     */
    public MCMoveGeometricClusterRestrictedGE(NeighborManager neighborManager1, NeighborManager neighborManager2,
                                              Space space, IRandom random, Box box1, Box box2, ISpecies seed, IPotential2[][] pairPotentials) {

        super();
        clusterAtoms1 = new HashSet<>();
        clusterAtoms2 = new HashSet<>();
        clusterAtomsList1 = new ArrayList<>();
        clusterAtomsList2 = new ArrayList<>();
        jNeighbors = new ArrayList<>();
        jNbrVecs = new ArrayList<>();
        pairNbrVecs = new ArrayList<>();
        this.neighborManager1 = neighborManager1;
        this.neighborManager2 = neighborManager2;
        this.neighborIterator1 = neighborManager1.makeNeighborIterator();
        this.neighborIterator2 = neighborManager2.makeNeighborIterator();
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
        potentials = pairPotentials;
        originalBox = new HashMap<>();
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
            for (int x = 0; x<jNeighbors.size(); x++) {
                IAtom atomJ = jNeighbors.get(x);
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                originalBox.put(atomJ, boxI);
                pairNbrVecs.add(jNbrVecs.get(x));
            }
            jNeighbors.clear();
            jNbrVecs.clear();
            moveAtom(atomI,boxI);
            Box otherBoxI = boxI==box1?box2:box1;
            gatherNeighbors(atomI, otherBoxI);
            for (int x = 0; x<jNeighbors.size(); x++) {
                IAtom atomJ = jNeighbors.get(x);
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                originalBox.put(atomJ, otherBoxI);
                pairNbrVecs.add(jNbrVecs.get(x));
            }
            while (atomPairs.size() > 0) {
                IAtom atomJ = atomPairs.remove(atomPairs.size() - 1);
                atomI = atomPairs.remove(atomPairs.size() - 1);
                Vector drij = pairNbrVecs.remove(pairNbrVecs.size() - 1);
                boxI = originalBox.get(atomI);
                Box boxJ = originalBox.get(atomJ);
                if(clusterAtoms2.contains(atomJ) || clusterAtoms1.contains(atomJ)) continue;
                double E;
                if(boxJ != boxI){
                    E = computeEnergy(drij, atomI, atomJ);
                }
                else{
                    moveAtom(atomJ, boxJ);
                    E = -1*computeEnergy(drij, atomI, atomJ);
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

    //Computes the energy between atomI and atomJ, determines the potential as needed and saves for later use.
    private double computeEnergy(Vector dr12, IAtom atomI, IAtom atomJ){
        int iIndex = atomI.getType().getIndex();
        int jIndex = atomJ.getType().getIndex();
        IPotential2 p = potentials[iIndex][jIndex];
        if(p == null) {
            return 0;
        }
        return p.u(dr12, atomI, atomJ);
    }

    private void gatherNeighbors(IAtom atomI, Box boxI) {
        if (boxI == box1) {
            neighborIterator1.iterAllNeighbors(atomI.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    if (clusterAtoms1.contains(jAtom)) return;
                    jNeighbors.add(jAtom);
                    Vector foo = boxI.getSpace().makeVector();
                    foo.E(rij);
                    jNbrVecs.add(foo);
                }
            });
        }
        else {
            neighborIterator2.iterAllNeighbors(atomI.getLeafIndex(), new NeighborIterator.NeighborConsumer() {
                @Override
                public void accept(IAtom jAtom, Vector rij) {
                    if (clusterAtoms2.contains(jAtom)) return;
                    jNeighbors.add(jAtom);
                    Vector foo = boxI.getSpace().makeVector();
                    foo.E(rij);
                    jNbrVecs.add(foo);
                }
            });
        }
    }

    public void setTemperature(double temperature){
        this.temperature = temperature;
    }

    private void moveAtom(IAtom atom, Box box){
        Vector position = atom.getPosition();
        position.TE(-1);
        position.PEa1Tv1(2,pivot);
        position.PE(box.getBoundary().centralImage(position));
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

    /**
     * Returns 1 because move is always accepted.
     */
    @Override
    public double energyChange(Box box) {
        return 0;
    }
}
