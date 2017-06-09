/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial;

import etomica.atom.IAtomList;
import etomica.atom.IMolecule;
import etomica.atom.IMoleculeList;
import etomica.atom.MoleculePositionGeometricCenter;
import etomica.space.Vector;
import etomica.atom.IMoleculePositionDefinition;
import etomica.space.Space;

public class ClusterCoupledFlippedPartial implements ClusterAbstract {

	public ClusterCoupledFlippedPartial(ClusterAbstract cluster, Space space, int[][] flipList) {
		this.space = space;
		this.flipList = flipList;
		actualFlipList = new int[flipList.length];
        wrappedCluster = cluster;
        childAtomVector = space.makeVector();
        flippedAtoms = new boolean[flipList.length];
        positionDefinition = new MoleculePositionGeometricCenter(space);
    }

	public ClusterAbstract makeCopy() {
		return new ClusterCoupledFlippedPartial(wrappedCluster.makeCopy(), space, flipList);
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

        for (int i=0; i<flipList.length; i++) {
            flippedAtoms[i] = false;
        }
        
        double vsum = wrappedCluster.value(box);
        int flipMol = 0;
        for (int i=0;i<flipList.length;i++){
        	int mol0 = flipList[i][0];
        	int mol1 = flipList[i][1];
        	
        	if (mol0 > mol1){
        		int swap = mol0;
        		mol0 = mol1;
        		mol1 = swap;
        	}
        	if(cPairs.getr2(mol0,mol1) > 100){//only flip molecules when the distance from the molecule0 is larger than 10A
        		actualFlipList[flipMol] = i;
        		flipMol++;
        	}
        }
        
        IMoleculeList atomList = box.getMoleculeList();
        // loop through the atoms, toggling each one until we toggle one "on"
        // this should generate each combination of flipped/unflipped for all
        // the molecules
        
        while (true) {
            boolean didFlipTrue = false;
            for (int i = 0; !didFlipTrue && i<flipMol; i++) {
                flippedAtoms[i] = !flippedAtoms[i];
                didFlipTrue = flippedAtoms[i];
                flip(atomList.getMolecule(flipList[actualFlipList[i]][0]));
                for (int j = 2;j < flipList[actualFlipList[i]].length; j++){
                	flip(atomList.getMolecule(flipList[actualFlipList[i]][j]),atomList.getMolecule(flipList[actualFlipList[i]][0]));
                }
            }
            if (!didFlipTrue) {
            	// if we flipped every atom from true to false, we must be done
                break;
            }
            vsum += wrappedCluster.value(box);
        }
        
        value = vsum / (1<<flipMol);
        
        cPairID = cPairs.getID();
        return value;
    }
    private void flip(IMolecule flippedMolecule, IMolecule centralMolecule) {
        Vector COM = positionDefinition.position(centralMolecule);
		IAtomList childAtoms = flippedMolecule.getChildList();
		for (int i = 0; i < childAtoms.getAtomCount(); i++) {
		    childAtomVector.Ea1Tv1(2,COM);
			childAtomVector.ME(childAtoms.getAtom(i).getPosition());
			childAtoms.getAtom(i).getPosition().E(childAtomVector);
		}
    }
    private void flip(IMolecule flippedMolecule) {
    	flip(flippedMolecule,flippedMolecule);
    }

    public void setTemperature(double temperature) {
        wrappedCluster.setTemperature(temperature);
    }
    
    private final ClusterAbstract wrappedCluster;
    protected final Space space;
    protected long cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected final boolean[] flippedAtoms;
    private Vector childAtomVector;
    protected IMoleculePositionDefinition positionDefinition;
    protected final int[][] flipList;
    protected final int[] actualFlipList;
}
