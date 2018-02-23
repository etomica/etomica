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

public class ClusterCoupledFlippedMultivalue implements ClusterAbstractMultivalue {

    /**
     * cluster must have caching disabled
     */
    public ClusterCoupledFlippedMultivalue(ClusterAbstractMultivalue cluster, Space space, int nDer) {
        this(cluster, space, 0,nDer);
    }
    
    /**
     * cluster must have caching disabled
     * configurations will be flipped when the minimum distance between any two molecules
     * exceeds minFlipDistance.  set minFlipDistance to 0 to always flip.
     */
    public ClusterCoupledFlippedMultivalue(ClusterAbstractMultivalue cluster, Space space, double minFlipDistance, int nDer) {
        this.space = space;
        wrappedCluster = cluster;
        childAtomVector = space.makeVector();
        flippedAtoms = new boolean[cluster.pointCount()];
        positionDefinition = new MoleculePositionGeometricCenter(space);
        this.minFlipDistance = minFlipDistance;
        value = new double[nDer+1];
        lastValue = new double[nDer+1];
        this.nDer = nDer;
        }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledFlippedMultivalue((ClusterAbstractMultivalue)wrappedCluster.makeCopy(), space, minFlipDistance, nDer);
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
        boolean debugme = false;
//        if (debugme)System.out.println("1 "+thisCPairID+" "+cPairID+" "+lastCPairID+" "+value[0]+" "+lastValue[0]+" "/*+f[0].getClass()*/);
        if (thisCPairID == cPairID) {
//            if (debugme)System.out.println("2 "+"clusterSum "+cPairID+" returning RECENT "+value[0]);
            return value[0];
        }
        else if (thisCPairID == lastCPairID) {
            // we went back to the previous cluster, presumably because the last
            // cluster was a trial that was rejected.  so drop the most recent value/ID
            cPairID = lastCPairID;
            for(int m=0;m<value.length;m++){
                value[m] = lastValue[m];
            }
//            if (debugme)System.out.println("3 "+"clusterSum "+cPairID+" returning PREVIOUS RECENT "+lastValue[0]);
            return value[0];
        }

        // a new cluster    
        lastCPairID = cPairID;      
        for(int m=0; m<value.length; m++){
            lastValue[m] = value[m];
        }
        cPairID = thisCPairID;

        final int pointCount = wrappedCluster.pointCount();
        
        boolean flipit = false;
        totcount++;
        double minR2 = minFlipDistance*minFlipDistance;
        double maxR2 = 0;
        for (int i=0; i<pointCount; i++) {
            flippedAtoms[i] = false;
            for (int j=i+1; j<pointCount; j++) {
                if (box.getCPairSet().getr2(i,j) > minR2) {
                    if (debugme && box.getCPairSet().getr2(i,j)>maxR2) maxR2 = box.getCPairSet().getr2(i,j);
                    if (false && box.getCPairSet().getr2(i,j) > 2*minR2) debugme=true; 
                    flipit=true; 
                    if(countflips)flipcount+=1;
                }
            }
        }
//        if (debugme)System.out.println("FLIPIT " + flipit);
        double[] v = wrappedCluster.getAllLastValues(box);
        for(int m=0;m<value.length;m++){
            value[m] = v[m];
        }
        
        if (!flipit) { 
//            if (debugme)System.out.println("4 "+"clusterSum "+cPairID+" returning V[0] "+v[0]);
            return v[0];
        }
        if (debugme) System.out.print(String.format("%6.3f %10.4e ", Math.sqrt(maxR2), v[0]));

        IMoleculeList atomList = box.getMoleculeList();
        // loop through the atoms, toggling each one until we toggle one "on"
        // this should generate each combination of flipped/unflipped for all
        // the molecules
        while (true) {
            boolean didFlipTrue = false;
            for (int i=0; !didFlipTrue && i<pointCount; i++) {
                flippedAtoms[i] = !flippedAtoms[i];
                didFlipTrue = flippedAtoms[i];
                flip(atomList.get(i));
            }
            if (!didFlipTrue) {
                // if we flipped every atom from true to false, we must be done
                break;
            }
            v = wrappedCluster.getAllLastValues(box);
            if (debugme) System.out.print(String.format("%10.4e ", v[0]));
            for(int m=0;m<value.length;m++){
                value[m] += v[m];
            }
            if (Double.isNaN(value[0])) {throw new RuntimeException("oops");}
        }
        for(int m=0;m<value.length;m++){
            value[m] /=  Math.pow(2, pointCount);
        }
        if (debugme) System.out.print(String.format("%10.4e\n", value[0]));
        
        
//        if (debugme) System.out.println("5 "+"clusterSum "+cPairID+" returning NEW "+value[0]);
        return value[0];
    }
    
    public double[] getAllLastValues(BoxCluster box) {
        value(box);
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
    
    public long getflipcount(){
        return flipcount;
    }
    
    public long gettotcount(){
        return totcount;
    }
    
    public double getflipfrac(){
        return ((double)flipcount)/totcount;
    }
    
    protected final ClusterAbstractMultivalue wrappedCluster;
    protected final Space space;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value[], lastValue[];
    protected final boolean[] flippedAtoms;
    private Vector childAtomVector;
    protected IMoleculePositionDefinition positionDefinition;
    protected final double minFlipDistance;
    protected final int nDer;
    protected final boolean countflips = true;
    protected long flipcount = 0;
    protected long totcount = 0;
}
