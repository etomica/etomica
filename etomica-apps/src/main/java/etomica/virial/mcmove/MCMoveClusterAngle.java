package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.data.meter.MeterPotentialEnergy;
import etomica.integrator.IntegratorMC;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.MoleculePositionCOM;
import etomica.molecule.MoleculePositionCOMPBC;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.potential.compute.PotentialComputeAggregate;
import etomica.space.RotationTensor;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.RotationTensor3D;
import etomica.space3d.Vector3D;
import etomica.species.ISpecies;
import etomica.util.collections.IntArrayList;
import etomica.util.random.IRandom;
import etomica.virial.BoxCluster;
import etomica.virial.PotentialComputeIntramolecular;

public class MCMoveClusterAngle extends MCMoveBoxStep {
    private final PotentialCompute potential;
    protected final IRandom random;
    protected final Space space;
    protected ISpecies species;
    protected Vector[][] position = null;
    protected final IntArrayList [] bonding;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    int [] modified;
    int modifiedIndex = 0;
    int b = 0;

    public MCMoveClusterAngle(PotentialCompute potentialCompute, Space space, IntArrayList[] bonding, IRandom random, double stepSize) {
        super(null);
        this.potential = potentialCompute;
        this.space = space;
        this.random = random;
        this.bonding = bonding;
        this.stepSize = stepSize;
        modified = new int[bonding.length];
    }

    public void setBox(Box p) {
        super.setBox(p);
        position = new Vector[p.getMoleculeList().size()][p.getMoleculeList().get(0).getChildList().size()];
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
        uOld = potential.computeAll(false);
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IMoleculeList moleculeList = box.getMoleculeList();
        for(int i = 0; i<moleculeList.size(); i++) {
            Vector com = space.makeVector();
            if (species != null && moleculeList.get(i).getType() != species) {
                continue;
            }
            IMolecule molecule = moleculeList.get(i);
            for(int j = 0; j < molecule.getChildList().size(); j++) {
                position[i][j] = space.makeVector();
                position[i][j].E(molecule.getChildList().getAtoms().get(j).getPosition());
            }
            com.E(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
            modifiedIndex = 0;
            double dt = stepSize * random.nextDouble();
            do{
                b = random.nextInt(bonding.length);
            }while (bonding[b].size() < 2);
            modified[modifiedIndex]=b;
            ++modifiedIndex;
            int a = random.nextInt(bonding[b].size());
            modified[modifiedIndex] = bonding[b].getInt(a);
            ++modifiedIndex;
            Vector axis = space.makeVector();
            axis.setRandomSphere(random);
            Vector r = space.makeVector();
            r.Ev1Mv2(molecule.getChildList().getAtoms().get(b).getPosition(), molecule.getChildList().getAtoms().get(bonding[b].getInt(a)).getPosition());
            Vector projection = space.makeVector();
            projection.Ea1Tv1(axis.dot(r)/r.squared(), r);
            axis.ME(projection);
            axis.normalize();
            RotationTensor3D rotationTensor = new RotationTensor3D();
            rotationTensor.setRotationAxis(axis, dt);
            transform(rotationTensor, bonding[b].getInt(a), molecule);
            transformBondedAtoms(rotationTensor, bonding[b].getInt(a), molecule);
            com.ME(MoleculePositionCOMPBC.com(box.getBoundary(), molecule));
            for(int j = 0; j < molecule.getChildList().size(); j++) {
                molecule.getChildList().getAtoms().get(j).getPosition().PE(com);
            }
        }
        ((BoxCluster)box).trialNotify();
        uNew = potential.computeAll(false);
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    protected void transformBondedAtoms(RotationTensor3D rotationTensor3D, int index, IMolecule molecule){
        for(int k = 0; k < bonding[index].size(); k++){
            boolean rotated = false;
            for (int l = 0; l < modifiedIndex; l++) {
                if (bonding[index].getInt(k) == modified[l]) {
                    rotated = true;
                    break;
                }
            }
            if (!rotated) {
                transform(rotationTensor3D, bonding[index].getInt(k), molecule);
                modified[modifiedIndex] = bonding[index].getInt(k);
                ++modifiedIndex;
                transformBondedAtoms(rotationTensor3D, bonding[index].getInt(k), molecule);
            }
        }
    }

    protected void transform(RotationTensor3D rotationTensor3D, int index, IMolecule molecule) {
        Vector r = space.makeVector();
        r.Ev1Mv2(molecule.getChildList().getAtoms().get(index).getPosition(), molecule.getChildList().getAtoms().get(b).getPosition());
        rotationTensor3D.transform(r);
        r.PE(molecule.getChildList().getAtoms().get(b).getPosition());
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
}
