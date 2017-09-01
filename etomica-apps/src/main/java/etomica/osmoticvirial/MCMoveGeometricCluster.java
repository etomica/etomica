package etomica.osmoticvirial;

import etomica.atom.*;
import etomica.atom.iterator.ApiBuilder;
import etomica.atom.iterator.AtomIterator;
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
import etomica.potential.Potential;
import etomica.potential.PotentialArray;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.Arrays;
import etomica.util.random.IRandom;

import java.util.HashSet;
import java.util.Random;

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

    public MCMoveGeometricCluster(PotentialMasterCell potentialMaster, Space space, IRandom random,
                                  double neighborRange, IntegratorMC integratorMC) {
        super(potentialMaster);
        oldPosition = space.makeVector();
        positionSource = new RandomPositionSourceRectangular(space, random);
        atomSource = new AtomSourceRandomLeaf();
        ((AtomSourceRandomLeaf)atomSource).setRandomNumberGenerator(random);
        clusterAtoms = new HashSet<>();
        jNeighbors = new HashSet<>();
        neighbors = new Api1ACell(space.D(),neighborRange,potentialMaster.getCellAgentManager());
        atomPairs = new AtomArrayList();
        neighbors.setDirection(null);
        this.integratorMC = integratorMC;
        this.random = random;
        potentials = new IPotentialAtomic[0][0];
        atomPair = new AtomPair();
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }

    @Override
    public boolean doTrial() {
        pivot = positionSource.randomPosition();
        IAtom atomI = atomSource.getAtom();
        if (atomI == null) return false;
        clusterAtoms.clear();
        NeighborCellManager ncm = ((PotentialMasterCell) potential).getNbrCellManager(box);
        outer: while(true){
            jNeighbors.clear();
            gatherNeighbors(atomI);
            moveAtom(atomI.getPosition());
            ncm.assignCell(atomI);
            gatherNeighbors(atomI);
            for (IAtom atomJ:jNeighbors) {
                atomPairs.add(atomI);
                atomPairs.add(atomJ);
            }
            while(atomPairs.getAtomCount() > 0){
                IAtom atomJ = atomPairs.remove(atomPairs.getAtomCount()-1);
                atomI = atomPairs.remove(atomPairs.getAtomCount()-1);
                double E = computeEnergy(atomI, atomJ);
                moveAtom(atomJ.getPosition());
                E -= computeEnergy(atomI, atomJ);
                double p = 1-Math.exp(-E/integratorMC.getTemperature());
                moveAtom(atomJ.getPosition());
                if(p<0 || p<random.nextDouble()) continue;
                atomI = atomJ;
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

    private void gatherNeighbors(IAtom atomI) {
        neighbors.setTarget(atomI);
        neighbors.reset();
        for(IAtomList pair = neighbors.next(); pair!= null; pair = neighbors.next()){
            IAtom atomJ = pair.getAtom(0);
            if(atomJ==atomI) atomJ = pair.getAtom(1);
            if (clusterAtoms.contains(atomJ)) continue;
            jNeighbors.add(atomJ);
        }
    }

    public void setBox(Box box){
        super.setBox(box);
        neighbors.setBox(box);
        positionSource.setBox(box);
    }

    private void moveAtom(Vector position){
        position.TE(-1);
        position.PEa1Tv1(2,pivot);
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
}
