/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.iterator.AtomIterator;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveBox;
import etomica.potential.PotentialMaster;

/**
 * Move that attempts to perform changes in the cluster diagram.
 */
public class MCMoveClusterDiagram extends MCMoveBox {

    public MCMoveClusterDiagram(PotentialMaster potentialMaster) {
        super(potentialMaster);
    }
    
    public void setBox(Box p) {
        super.setBox(p);
        ClusterAbstract sampleCluster = ((BoxCluster)p).getSampleCluster();
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
    }
    
    public boolean doTrial() {
        // don't notify the box.  we're not moving any atoms.
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        cluster.randomizeDiagram();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    public double getChi(double temperature) {
//        System.out.println("uNew "+uNew+" uOld "+uOld);
        foo += uNew * uOld;
        foo2 += uOld * uOld;
        return uNew/uOld;
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
