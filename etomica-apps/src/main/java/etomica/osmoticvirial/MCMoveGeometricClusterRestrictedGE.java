package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.atom.iterator.AtomIterator;
import etomica.atom.iterator.AtomIteratorArrayListSimple;
import etomica.box.Box;
import etomica.box.RandomPositionSource;
import etomica.box.RandomPositionSourceRectangular;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMove;
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
import etomica.util.Arrays;
import etomica.util.random.IRandom;

import java.util.ArrayList;
import java.util.HashSet;

public class MCMoveGeometricClusterRestrictedGE extends MCMove {

    protected RandomPositionSource positionSource;
    protected AtomSource atomSource;
    protected final HashSet<IAtom> clusterAtoms1, clusterAtoms2, jNeighbors;
    protected final Vector oldPosition;
    protected Vector pivot;
    protected final Api1ACell neighbors;
    protected final AtomArrayList atomPairs;
    protected final ArrayList<Box> boxArrayList;
    protected final IRandom random;
    protected IPotentialAtomic[][] potentials;
    protected final AtomPair atomPair;
    protected final AtomIteratorArrayListSimple atomIterator;
    protected final Box box1, box2;
    protected double temperature;

    public MCMoveGeometricClusterRestrictedGE(PotentialMasterCell potentialMaster, Space space, IRandom random,
                                              double neighborRange, Box box1, Box box2) {

        super(potentialMaster);
        clusterAtoms1 = new HashSet<>();
        clusterAtoms2 = new HashSet<>();
        jNeighbors = new HashSet<>();
        neighbors = new Api1ACell(space.D(),neighborRange,potentialMaster.getCellAgentManager());
        atomPairs = new AtomArrayList();
        oldPosition = space.makeVector();
        positionSource = new RandomPositionSourceRectangular(space, random);
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        neighbors.setDirection(null);
        this.box1 = box1;
        this.box2 = box2;
        this.random = random;
        potentials = new IPotentialAtomic[0][0];
        atomPair = new AtomPair();
        atomIterator = new AtomIteratorArrayListSimple();
        boxArrayList = new ArrayList<>();
    }

    @Override
    public boolean doTrial() {
        positionSource.setBox(box1);
        pivot = positionSource.randomPosition();
        atomSource.setBox(box1);
        Box boxI = box1; //change this to random in future
        IAtom atomI = atomSource.getAtom();
        if (atomI == null) return false;
        System.out.println("in geo cluster" + moveTracker);
        System.out.println(atomI + "initial"+boxI);
        clusterAtoms1.clear();
        clusterAtoms2.clear();
        clusterAtoms2.add(atomI);
       // System.out.println("molecule count "+box.getMoleculeList().getMoleculeCount());
        outer: while(true){
            jNeighbors.clear();
            gatherNeighbors(atomI, boxI);
            for (IAtom atomJ:jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                boxArrayList.add(boxI);
            }
            jNeighbors.clear();
            moveAtom(atomI,boxI);
            Box otherBoxI = boxI==box1?box2:box1;
            gatherNeighbors(atomI, otherBoxI);
            for (IAtom atomJ:jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
                boxArrayList.add(otherBoxI);
                 }
            while(atomPairs.getAtomCount() > 0){
                IAtom atomJ = atomPairs.remove(atomPairs.getAtomCount()-1);
                atomI = atomPairs.remove(atomPairs.getAtomCount()-1);
                Box boxJ = boxArrayList.remove(boxArrayList.size()-1);
                if(clusterAtoms2.contains(atomJ) || clusterAtoms1.contains(atomJ)) continue;
                double E = 0;
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
//                System.out.println(atomI + "initial"+boxI);
                HashSet<IAtom> clusterAtoms = boxI==box1?clusterAtoms2:clusterAtoms1;
                clusterAtoms.add(atomI);
                continue outer;
            }
            break;
        }

        return true;

    }

    private double computeEnergy(IAtom atomI, IAtom atomJ){
        int iIndex = atomI.getType().getIndex();
        int jIndex = atomJ.getType().getIndex();
        if(potentials.length <= iIndex){
            int oldSize = potentials.length;
            potentials = (IPotentialAtomic[][]) Arrays.resizeArray(potentials, iIndex+1);
            for(int i=oldSize; i<potentials.length; i++){
                potentials[i] = new IPotentialAtomic[0];
            }
        }
        if(potentials[iIndex].length <= jIndex) potentials[iIndex] = (IPotentialAtomic[]) Arrays.resizeArray(potentials[iIndex], jIndex+1);
        IPotentialAtomic p = potentials[iIndex][jIndex];
        atomPair.atom0 = atomI;
        atomPair.atom1 = atomJ;
        if(p == null) {
            PotentialArray iPotentials = ((PotentialMasterCell)potential).getRangedPotentials(atomI.getType());
            NeighborCriterion[] neighborCriteria = iPotentials.getCriteria();
            IPotential[] iPotential = iPotentials.getPotentials();
            for(int i=0; i<iPotential.length; i++){
                if(neighborCriteria[i].accept(atomPair)) {
                    potentials[iIndex][jIndex] = (IPotentialAtomic) iPotential[i];
                    p = potentials[iIndex][jIndex];
                }
            }
        }
        return p.energy(atomPair);
    }

    private void gatherNeighbors(IAtom atomI, Box boxI) {
        neighbors.setTarget(atomI);
        neighbors.setBox(boxI);
        neighbors.reset();
        for(IAtomList pair = neighbors.next(); pair!= null; pair = neighbors.next()){
            IAtom atomJ = pair.getAtom(0);
            if(atomJ==atomI) atomJ = pair.getAtom(1);
            if (clusterAtoms1.contains(atomJ) || clusterAtoms2.contains(atomJ)) continue;
            jNeighbors.add(atomJ);
        }
    }


    public void setTemperature(double temperature){
        this.temperature = temperature;
    }

    private void moveAtom(IAtom atom, Box box){
        Vector position = atom.getPosition();
        position.TE(-1);
        position.PEa1Tv1(2,pivot);
//        System.out.println(" remove from "+box+" "+atom+" "+atom.hashCode());
        box.removeMolecule(atom.getParentGroup());
        Box otherBox = box==box1?box2:box1;
        otherBox.addMolecule(atom.getParentGroup());
//        System.out.println(" add to "+otherBox+" "+atom+" "+atom.hashCode());
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
        HashSet<IAtom> clusterAtoms = box==box1?clusterAtoms1:clusterAtoms2;
        for(IAtom atom: clusterAtoms){
            atomList.add(atom);
        }
        return atomIterator;
    }

    @Override
    public double energyChange(Box box) {
        return 0;
    }
}
