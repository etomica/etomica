/*
 * Created on May 16, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.Phase;
import etomica.PotentialMaster;
import etomica.Simulation;
import etomica.Space;
import etomica.action.AtomActionTransform;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MCMoveClusterRotateMolecule3D extends MCMoveRotateMolecule3D implements MCMoveCluster {

    private final MeterClusterWeight weightMeter;

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

    public boolean doTrial() {
        Phase phase = phases[0];
        if(phase.moleculeCount()==1) {molecule = null; return false;}
            
        molecule = phase.randomMolecule();
        while (molecule.node.getOrdinal() == 1) {
            molecule = phase.randomMolecule();
        }
        uOld = weightMeter.getDataAsScalar(phase);
        
        double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(dTheta);

        leafAtomIterator.setRoot(molecule);
        leafAtomIterator.reset();
        r0.E(molecule.node.firstLeafAtom().coord.position());
        AtomActionTransform.doAction(leafAtomIterator, r0, rotationTensor);
            
        uNew = Double.NaN;
        ((PhaseCluster)phases[0]).trialNotify();
        return true;
    }
    
    public double trialRatio() {return 1.0;}
    
    public double lnProbabilityRatio() {
        return Math.log(probabilityRatio());
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar(phases[0]);
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
        
}
