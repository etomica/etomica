/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.atom.IAtomList;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculePositionDefinition;
import etomica.molecule.MoleculePositionGeometricCenter;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.virial.BoxCluster;
import etomica.virial.CoordinatePairSet;

import java.util.ArrayList;

public class ClusterCoupledFlippedPoints implements ClusterAbstract {

    private final int[][] flipPoints;

    /**
     * cluster must have caching disabled
     */
    public ClusterCoupledFlippedPoints(ClusterAbstract cluster, Space space, int[][] flipPoints) {
        this(cluster, space, flipPoints, 0);
    }
    
    /**
     * cluster must have caching disabled
     * configurations will be flipped when the minimum distance between any two molecules
     * exceeds minFlipDistance.  set minFlipDistance to 0 to always flip.
     */
    public ClusterCoupledFlippedPoints(ClusterAbstract cluster, Space space, int[][] flipPoints, double minFlipDistance) {
        this.space = space;
        wrappedCluster = cluster;
        childAtomVector = space.makeVector();
        flippedAtoms = new boolean[cluster.pointCount()];
        positionDefinition = new MoleculePositionGeometricCenter(space);
        this.minFlipDistance = minFlipDistance;
        this.flipPoints = flipPoints;
    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledFlippedPoints(wrappedCluster.makeCopy(), space, flipPoints, minFlipDistance);
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
        
        double minR2 = minFlipDistance*minFlipDistance;
        boolean debugme = false;
        ArrayList<int[]> myflipPoints = new ArrayList<>();
        for (int i=0; i<flipPoints.length; i++) {
            flippedAtoms[i] = false;
            if (box.getCPairSet().getr2(flipPoints[i][0],flipPoints[i][1]) > minR2) {
                myflipPoints.add(flipPoints[i]);
            }
        }
        
        double vsum = wrappedCluster.value(box);
        if (myflipPoints.size() == 0) {
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
            for (int i=0; !didFlipTrue && i<myflipPoints.size(); i++) {
                flippedAtoms[i] = !flippedAtoms[i];
                didFlipTrue = flippedAtoms[i];
                flip(atomList, myflipPoints.get(i));
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
        
        value = vsum / Math.pow(2, myflipPoints.size());
        if (debugme) System.out.print(String.format("%10.4e\n", value));
        
        return value;
    }
    
    protected void flip(IMoleculeList molecules, int[] myflipPoints) {
        Vector COM = positionDefinition.position(molecules.get(myflipPoints[1]));
        for (int j=1; j<myflipPoints.length; j++) {
            IAtomList childAtoms = molecules.get(j).getChildList();
            for (int i = 0; i < childAtoms.size(); i++) {
                childAtomVector.Ea1Tv1(2,COM);
                childAtomVector.ME(childAtoms.get(i).getPosition());
                childAtoms.get(i).getPosition().E(childAtomVector);
            }

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
