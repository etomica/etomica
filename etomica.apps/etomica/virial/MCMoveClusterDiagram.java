/*
 * Created on Oct 1, 2005
 *
 * TODO To change the template for this generated file go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
package etomica.virial;

import etomica.atom.iterator.AtomIterator;
import etomica.integrator.mcmove.MCMovePhase;
import etomica.phase.Phase;
import etomica.potential.PotentialMaster;

/**
 * @author andrew
 *
 * TODO To change the template for this generated type comment go to
 * Window - Preferences - Java - Code Style - Code Templates
 */
public class MCMoveClusterDiagram extends MCMovePhase {

    private MeterClusterWeight weightMeter;
    public MCMoveClusterDiagram(PotentialMaster potentialMaster) {
        super(potentialMaster);
        weightMeter = new MeterClusterWeight(potentialMaster);
    }
    
    public void setPhase(Phase p) {
        super.setPhase(p);
        weightMeter.setPhase(p);
        ClusterAbstract sampleCluster = ((PhaseCluster)p).getSampleCluster();
        ClusterAbstract cluster1;
        if (sampleCluster instanceof ClusterWeightAbs) {
            cluster1 = ((ClusterWeightAbs)sampleCluster).getSubCluster();
        }
        else { // must be umbrella
            cluster1 = ((ClusterWeightUmbrella)sampleCluster).getClusters()[0];
        }
        if (cluster1 instanceof ClusterSumStickyEF) {
            cluster = (ClusterSumStickyEF)cluster1;
        }
        else {
            cluster = (ClusterSumStickyEF)((ClusterCoupled)cluster1).getSubCluster();
        }
    }
    
    public boolean doTrial() {
        // don't notify the phase.  we're not moving any atoms.
        uOld = weightMeter.getDataAsScalar();
        cluster.randomizeDiagram();
        return true;
    }
    
    public double getA() {
        uNew = weightMeter.getDataAsScalar();
//        System.out.println("uNew "+uNew+" uOld "+uOld);
        foo += uNew * uOld;
        foo2 += uOld * uOld;
        return uNew/uOld;
    }
    
    public double getB() {
        return 0;
    }
    
    public void acceptNotify() {
  //      System.out.println("accepted"+(uNew == uOld ? " no change" : ""));
        // do nothing
    }
    
    public void rejectNotify() {
//        System.out.println("rejected");
        cluster.revertDiagram();
    }
    
    public double energyChange() {
        return uNew/uOld;
    }
    
    public AtomIterator affectedAtoms() {
        // you deserve it.
        return null;
    }
    
    private double uOld, uNew;
    private ClusterSumStickyEF cluster;
    public static double foo, foo2;
}
