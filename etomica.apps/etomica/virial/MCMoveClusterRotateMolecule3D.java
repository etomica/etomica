package etomica.virial;

import etomica.action.AtomAction;
import etomica.action.AtomTransform;
import etomica.atom.AtomLeaf;
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
public class MCMoveClusterRotateMolecule3D extends MCMoveRotateMolecule3D {

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
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
        oldPositions = new Vector[((AtomTreeNodeGroup)molecule.getNode()).getChildList().size()-1];
        for (int j=0; j<oldPositions.length; j++) {
            oldPositions[j] = p.space().makeVector();
        }
    }

    public boolean doTrial() {
        if(phase.moleculeCount()==1) {molecule = null; return false;}
            
        molecule = phase.randomMolecule();
        while (molecule.getNode().getIndex() == 0) {
            molecule = phase.randomMolecule();
        }
        uOld = weightMeter.getDataAsScalar();
        
        double dTheta = (2*Simulation.random.nextDouble() - 1.0)*stepSize;
        rotationTensor.setAxial(Simulation.random.nextInt(3),dTheta);

        leafAtomIterator.setRoot(molecule);
        leafAtomIterator.reset();
        AtomLeaf first = (AtomLeaf)leafAtomIterator.nextAtom();
        int j=0;
        while (leafAtomIterator.hasNext()) {
            oldPositions[j++].E(((AtomLeaf)leafAtomIterator.nextAtom()).getCoord().position());
        }
        leafAtomIterator.reset();
        r0.E(first.getCoord().position());
        AtomTransform.doTransform(leafAtomIterator, r0, rotationTensor);
            
        if (trialCount-- == 0) {
            relaxAction.setAtom(molecule);
            relaxAction.actionPerformed();
            trialCount = relaxInterval;
        }

        uNew = Double.NaN;
        ((PhaseCluster)phase).trialNotify();
        return true;
    }
    
    public double getB() {
        return 0.0;
    }
    
    public double getA() {
        uNew = weightMeter.getDataAsScalar();
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }
    
    public void acceptNotify() {
        super.acceptNotify();
        ((PhaseCluster)phase).acceptNotify();
    }
    
    public void rejectNotify() {
        super.rejectNotify();
        ((PhaseCluster)phase).rejectNotify();
    }
    
    public void setRelaxAction(AtomAction action) {
        relaxAction = action;
    }
    
    private static final long serialVersionUID = 1L;
    private final MeterClusterWeight weightMeter;
    protected int trialCount, relaxInterval = 100;
    protected AtomAction relaxAction;
    private Vector[] oldPositions;

}
