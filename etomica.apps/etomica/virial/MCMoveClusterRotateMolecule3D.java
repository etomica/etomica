package etomica.virial;

import etomica.action.AtomTransform;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */

/*
 * Created on May 16, 2005
 */

public class MCMoveClusterRotateMolecule3D extends MCMoveRotateMolecule3D implements MCMoveCluster {

    /**
     * @param potentialMaster
     * @param space
     */
    public MCMoveClusterRotateMolecule3D(PotentialMaster potentialMaster,
            Space space) {
        super(potentialMaster, space);
        weightMeter = new MeterClusterWeight(potential);
        setName("MCMoveClusterMolecule");
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
    }

    public boolean doTrial() {
        Phase phase = phases[0];
        if(phase.moleculeCount()==1) {molecule = null; return false;}
            
        molecule = phase.randomMolecule();
        while (molecule.node.getOrdinal() == 1) {
            molecule = phase.randomMolecule();
        }
        uOld = weightMeter.getDataAsScalar();
        
        double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(Simulation.random.nextInt(3),dTheta);

        leafAtomIterator.setRoot(molecule);
        leafAtomIterator.reset();
        r0.E(molecule.node.firstLeafAtom().coord.position());
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
            
        uNew = Double.NaN;
        ((PhaseCluster)phases[0]).trialNotify(null);
        return true;
    }
    
    public double trialRatio() {return 1.0;}
    
    public double lnProbabilityRatio() {
        return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar();
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phases[0]).acceptNotify();
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        ((PhaseCluster)phases[0]).rejectNotify();
    }
        
    private final MeterClusterWeight weightMeter;

}
