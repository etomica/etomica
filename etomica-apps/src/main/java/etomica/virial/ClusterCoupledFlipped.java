/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Space;
import etomica.space.Vector;

public class ClusterCoupledFlipped implements ClusterAbstract {

    /**
     * cluster must have caching disabled
     */
    public ClusterCoupledFlipped(ClusterAbstract cluster, Space space) {
        this(cluster, space, 0);
    }
    
    /**
     * cluster must have caching disabled
     * configurations will be flipped when the minimum distance between any two molecules
     * exceeds minFlipDistance.  set minFlipDistance to 0 to always flip.
     */
    public ClusterCoupledFlipped(ClusterAbstract cluster, Space space, double minFlipDistance) {
        this.space = space;
        wrappedCluster = cluster;
        childAtomVector = space.makeVector();
        flippedAtoms = new boolean[cluster.pointCount()];
        positionDefinition = new MoleculePositionGeometricCenter(space);
        this.minFlipDistance = minFlipDistance;
    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledFlipped(wrappedCluster.makeCopy(), space, minFlipDistance);
    }

    public int pointCount() {
        return wrappedCluster.pointCount();
    }

    public ClusterAbstract getSubCluster() {
        return wrappedCluster;
    }
    
    public double value(BoxCluster box) {
        CoordinatePairSet cPairs = box.getCPairSet();
        long thisCPairID = cPairs.getID();
//      System.out.println(thisCPairID+" "+cPairID+" "+lastCPairID+" "+value+" "+lastValue+" "+f[0].getClass());
        if (thisCPairID == cPairID) {
//          System.out.println("clusterSum "+cPairID+" returning recent "+value);
            return value;
        }
        else if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            value = lastValue;
//          System.out.println("clusterSum "+cPairID+" returning previous recent "+lastValue);
            return value;
        }

        // a new cluster
        lastCPairID = cPairID;
        lastValue = value;
        cPairID = thisCPairID;

        final int pointCount = wrappedCluster.pointCount();
        
        boolean flipit = false;
        double minR2 = minFlipDistance*minFlipDistance;
        boolean debugme = false;
        for (int i=0; i<pointCount; i++) {
            flippedAtoms[i] = false;
            for (int j=i+1; j<pointCount; j++) {
                if (box.getCPairSet().getr2(i,j) > minR2) {
                    if (false && box.getCPairSet().getr2(i,j) > 2*minR2) debugme=true; 
                    flipit=true;
                }
            }
        }
        
        double vsum = wrappedCluster.value(box);
        if (!flipit) {
            value = vsum;
            return vsum;
        }
        if (debugme) System.out.print(String.format("%10.4e ", vsum));

        IMoleculeList atomList = box.getMoleculeList();
        // loop through the atoms, toggling each one until we toggle one "on"
        // this should generate each combination of flipped/unflipped for all
        // the molecules
        while (true) {
            boolean didFlipTrue = false;
            for (int i=0; !didFlipTrue && i<pointCount; i++) {
                flippedAtoms[i] = !flippedAtoms[i];
                didFlipTrue = flippedAtoms[i];
                flip(atomList.getMolecule(i));
            }
            if (!didFlipTrue) {
                // if we flipped every atom from true to false, we must be done
                break;
            }
            double foo = wrappedCluster.value(box);
            if (debugme) System.out.print(String.format("%10.4e ", foo));
            vsum += foo;
            if (Double.isNaN(vsum)) {throw new RuntimeException("oops");}
        }
        
        value = vsum / Math.pow(2, pointCount);
        if (debugme) System.out.print(String.format("%10.4e\n", value));
        
        return value;
    }
    
    protected void flip(IMolecule flippedMolecule) {
        Vector COM = positionDefinition.position(flippedMolecule);
		IAtomList childAtoms = flippedMolecule.getChildList();
		for (int i = 0; i < childAtoms.size(); i++) {
		    childAtomVector.Ea1Tv1(2,COM);
			childAtomVector.ME(childAtoms.get(i).getPosition());
			childAtoms.get(i).getPosition().E(childAtomVector);
		}
    }

    public void setTemperature(double temperature) {
        wrappedCluster.setTemperature(temperature);
    }
    
    protected final ClusterAbstract wrappedCluster;
    protected final Space space;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected final boolean[] flippedAtoms;
    private Vector childAtomVector;
    protected IMoleculePositionDefinition positionDefinition;
    protected final double minFlipDistance;
}
