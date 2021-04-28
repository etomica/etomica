package etomica.virial.mcmove;

import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.integrator.mcmove.MCMoveStep;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;

public class MCMoveClusterTorsion extends MCMoveBoxStep {

    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected ISpecies species;
    protected double dt = 0;
    protected int a,b;
    protected Vector[][] position = null;
    protected final IntArrayList [] bonding;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    int [] indices = new int[3];
    int [] modified;
    int modifiedIndex = 0;
    //int com = 0;

    public MCMoveClusterTorsion(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize) {
        super(null);
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.bonding = bonding;
        this.stepSize = stepSize;
        setStepSizeMin(stepSize);
        modified = new int[bonding.length];
    }

    public void setBox(Box p) {
        super.setBox(p);
        position = new Vector[p.getMoleculeList().size()][p.getMoleculeList().get(0).getChildList().size()];
    }

    @Override
    public boolean doTrial() {
        uOld = potential.computeAll(false);
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<moleculeList.size(); i++) {
            if (species != null && moleculeList.get(i).getType() != species) {
                continue;
            }
            IMolecule molecule = moleculeList.get(i);
            for(int j = 0; j < molecule.getChildList().size(); j++) {
                position[i][j] = space.makeVector();
                position[i][j].E(molecule.getChildList().getAtoms().get(j).getPosition());
            }
            modifiedIndex = 0;
            dt = stepSize * (random.nextDouble() - 0.5);
            do{
                a = random.nextInt(bonding.length);
            }while (bonding[a].size() < 2);
            modified[modifiedIndex]=a;
            ++modifiedIndex;
            do{
                b = bonding[a].getInt(random.nextInt(bonding[a].size()));
            }while(bonding[b].size() < 2);
            modified[modifiedIndex] = b;
            ++modifiedIndex;
            Vector axis = space.makeVector();
            axis.Ev1Mv2(molecule.getChildList().getAtoms().get(a).getPosition(), molecule.getChildList().getAtoms().get(b).getPosition());
            axis.normalize();
            RotationTensor3D rotationTensor = new RotationTensor3D();
            rotationTensor.setRotationAxis(axis, dt/2);
            transformBondedAtoms(rotationTensor, position[i][a], a, molecule);
            rotationTensor.invert();
            transformBondedAtoms(rotationTensor, position[i][a], b, molecule);
        }
        ((BoxCluster)box).trialNotify();
        uNew = potential.computeAll(false);
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    protected void transformBondedAtoms(RotationTensor3D rotationTensor3D, Vector pointOnAxis, int index, IMolecule molecule){
        for(int k = 0; k < bonding[index].size(); k++){
            boolean rotated = false;
            for (int l = 0; l < modifiedIndex; l++) {
                if (bonding[index].getInt(k) == modified[l]) {
                    rotated = true;
                    break;
                }
            }
            if (!rotated) {
                transform(rotationTensor3D, pointOnAxis, bonding[index].getInt(k), molecule);
                modified[modifiedIndex] = bonding[index].getInt(k);
                ++modifiedIndex;
                transformBondedAtoms(rotationTensor3D, pointOnAxis, bonding[index].getInt(k), molecule);
            }
        }
    }

    protected void transform(RotationTensor3D rotationTensor3D, Vector pointOnAxis, int index, IMolecule molecule) {
        //System.out.println("index " + index + " pointOnAxis " + pointOnAxis);
        Vector r = space.makeVector();
        r.Ev1Mv2(molecule.getChildList().getAtoms().get(index).getPosition(), pointOnAxis);
        rotationTensor3D.transform(r);
        r.PE(pointOnAxis);
        molecule.getChildList().getAtoms().get(index).getPosition().E(r);
    }

    @Override
    public double getChi(double temperature) {
        return (wOld == 0 ? 1 : wNew / wOld) * Math.exp(-(uNew - uOld) / temperature);
    }

    @Override
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    @Override
    public void rejectNotify() {
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<box.getMoleculeList().size(); i++) {
            for(int j = 0; j < moleculeList.get(i).getChildList().size(); j++) {
                moleculeList.get(i).getChildList().getAtoms().get(j).getPosition().E(position[i][j]);
            }
        }
        ((BoxCluster)box).rejectNotify();
    }

    @Override
    public AtomIterator affectedAtoms() {
        return null;
    }

    @Override
    public double energyChange() {
        return 0;
    }
}
