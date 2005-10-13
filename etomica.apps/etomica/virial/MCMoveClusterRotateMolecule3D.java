package etomica.virial;

import etomica.action.AtomAction;
import etomica.action.AtomTransform;
import etomica.atom.Atom;
import etomica.atom.AtomTreeNodeGroup;
import etomica.integrator.mcmove.MCMoveRotateMolecule3D;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;

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
        oldPositions = new Vector[((AtomTreeNodeGroup)molecule.node).childList.size()-1];
        for (int j=0; j<oldPositions.length; j++) {
            oldPositions[j] = p[0].space().makeVector();
        }
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
        Atom first = leafAtomIterator.nextAtom();
        int j=0;
        while (leafAtomIterator.hasNext()) {
            oldPositions[j++].E(leafAtomIterator.nextAtom().coord.position());
        }
        leafAtomIterator.reset();
        r0.E(first.coord.position());
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
            
        if (trialCount-- == 0) {
            relaxAction.setAtom(molecule);
            relaxAction.actionPerformed();
            trialCount = relaxInterval;
        }

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
    
    public void setRelaxAction(AtomAction action) {
        relaxAction = action;
    }
    
    private final MeterClusterWeight weightMeter;
    protected int trialCount, relaxInterval = 100;
    protected AtomAction relaxAction;
    private Vector[] oldPositions;

}
