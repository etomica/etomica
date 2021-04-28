package etomica.virial.mcmove;

import etomica.atom.IAtomList;
import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.IntegratorMDFasterer;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.integrator.mcmove.MCMoveBoxStep;
import etomica.potential.PotentialMaster;
import etomica.potential.compute.PotentialCompute;
import etomica.space.Vector;
import etomica.virial.BoxCluster;
import etomica.virial.PotentialComputeIntramolecular;

import java.util.List;

public class MCMoveClusterMD extends MCMoveBoxStep {

    IntegratorMDFasterer integratorMDFaster;
    double uOld = 0;
    double uNew = 0;
    double wOld = 0;
    double wNew = 0;
    Vector[] position = null;
    long steps = 20;
    public MCMoveClusterMD(IntegratorMDFasterer integratorMDFasterer, Box box) {
        super(null);
        setBox(box);
        this.integratorMDFaster = integratorMDFasterer;
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
        //PotentialComputeIntramolecular.doDebug = true;
        integratorMDFaster.reset();
        integratorMDFaster.doThermostat();
        //System.out.printf("Initial PE %f , KE %f, TE %f \n",integratorMDFaster.getPotentialEnergy(), integratorMDFaster.getKineticEnergy(), integratorMDFaster.getPotentialEnergy() + integratorMDFaster.getKineticEnergy());
        PotentialComputeIntramolecular.doDebug = false;
        IAtomList atomList = integratorMDFaster.getBox().getLeafList();
        position = new Vector[atomList.size()];
        //System.out.printf("Starting position");
        for (int i = 0; i < atomList.size(); i++) {
            position[i] = integratorMDFaster.getBox().getSpace().makeVector();
            position[i].E(atomList.get(i).getPosition());
            //System.out.print(position[i] + " ");
        }
        //System.out.println();
        uOld = integratorMDFaster.getPotentialEnergy() + integratorMDFaster.getKineticEnergy();
        wOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        for (int i = 0; i < steps; i++) {
            integratorMDFaster.doStep();
            //System.out.printf("Steps %d, PE %f , KE %f, TE %f \n", i, integratorMDFaster.getPotentialEnergy(), integratorMDFaster.getKineticEnergy(), integratorMDFaster.getPotentialEnergy() + integratorMDFaster.getKineticEnergy());
        }
        //System.exit(0);
        ((BoxCluster)box).trialNotify();
        uNew = integratorMDFaster.getPotentialEnergy() + integratorMDFaster.getKineticEnergy();
        wNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    @Override
    public double getChi(double temperature) {
        return Math.exp(-(uNew - uOld)/temperature)*(wNew/wOld);
    }

    @Override
    public void acceptNotify() {
        ((BoxCluster)box).acceptNotify();
    }

    @Override
    public void rejectNotify() {
        IAtomList atomList = integratorMDFaster.getBox().getLeafList();
        //System.out.printf("Revert to original position");
        for (int i = 0; i < atomList.size(); i++) {
            atomList.get(i).getPosition().E(position[i]);
            //System.out.print(position[i] + " ");
        }
        //System.out.println();

        ((BoxCluster)box).rejectNotify();
    }
}
