/*
 * Created on Oct 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.atom.iterator.AtomIterator;
import etomica.integrator.MCMove;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MCMoveClusterDiagram extends MCMove implements MCMoveCluster {

    private MeterClusterWeight weightMeter;
    public MCMoveClusterDiagram(PotentialMaster potentialMaster) {
        super(potentialMaster,1);
        weightMeter = new MeterClusterWeight(potentialMaster);
    }
    
    public void setPhase(Phase[] p) {
        super.setPhase(p);
        weightMeter.setPhase(p[0]);
        cluster = (ClusterSumSticky)((ClusterWeightAbs)((PhaseCluster)p[0]).getSampleCluster()).getSubCluster();
    }
    
    public boolean doTrial() {
        // don't notify the phase.  we're not moving any atoms.
        uOld = weightMeter.getDataAsScalar();
        cluster.randomizeDiagram();
        return true;
    }
    
    public double trialRatio() {
        return 1;
    }
    
    public double probabilityRatio() {
        uNew = weightMeter.getDataAsScalar();
//        System.out.println("uNew "+uNew+" uOld "+uOld);
        return uNew/uOld;
    }
    
    public double lnTrialRatio() {
        return Math.log(trialRatio());
    }
    
    public double lnProbabilityRatio() {
        return 0;
    }
    
    public void acceptNotify() {
//        System.out.println("accepted");
        // do nothing
    }
    
    public void rejectNotify() {
//        System.out.println("rejected");
        cluster.revertDiagram();
    }
    
    public double energyChange(Phase phase) {
        return uNew/uOld;
    }
    
    public AtomIterator affectedAtoms(Phase phase) {
        // you deserve it.
        return null;
    }
    
    private double uOld, uNew;
    private ClusterSumSticky cluster;
}
