package etomica.virial;

import etomica.atom.AtomGroup;
import etomica.atom.AtomSet;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;
import etomica.space.IVector;
import etomica.space3d.Vector3D;

public class ClusterCoupledFlipped implements ClusterAbstract {

    public ClusterCoupledFlipped(ClusterAbstract cluster) {
        wrappedCluster = cluster;
        childAtomVector = new Vector3D();
        flippedAtoms = new boolean[cluster.pointCount()];
    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupledFlipped(wrappedCluster);
    }

    public int pointCount() {
        return wrappedCluster.pointCount();
    }

    public ClusterAbstract getSubCluster() {
        return wrappedCluster;
    }
    
    public double value(BoxCluster box) {
        AtomSet atomList = ((AtomGroup)box.getSpeciesMaster().getAgentList().getAtom(0)).getChildList();
        CoordinatePairSet cPairs = box.getCPairSet();
        int thisCPairID = cPairs.getID();
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
        
        for (int i=0; i<pointCount; i++) {
            flippedAtoms[i] = false;
        }
        
        double vsum = wrappedCluster.value(box);

        // loop through the atoms, toggling each one until we toggle one "on"
        // this should generate each combination of flipped/unflipped for all
        // the molecules
        while (true) {
            boolean didFlipTrue = false;
            for (int i=0; !didFlipTrue && i<pointCount; i++) {
                flippedAtoms[i] = !flippedAtoms[i];
                didFlipTrue = flippedAtoms[i];
                flip((IAtomGroup)atomList.getAtom(i));
            }
            if (!didFlipTrue) {
                // if we flipped every atom from true to false, we must be done
                break;
            }
            cPairs.reset();
            vsum += wrappedCluster.value(box);
        }
        
        value = vsum / Math.pow(2, pointCount);
        
        cPairID = cPairs.getID();
        return value;
    }
    
    private void flip(IAtomGroup flippedMolecule) {
        IVector COM = flippedMolecule.getType().getPositionDefinition().position(flippedMolecule);
		AtomSet childAtoms = flippedMolecule.getChildList();
		for (int i = 0; i < childAtoms.getAtomCount(); i++) {
		    childAtomVector.Ea1Tv1(2,COM);
			childAtomVector.ME(((IAtomPositioned)childAtoms.getAtom(i)).getPosition());
			((IAtomPositioned)childAtoms.getAtom(i)).getPosition().E(childAtomVector);
		}
    }

    public void setTemperature(double temperature) {
        wrappedCluster.setTemperature(temperature);
    }
    
    private final ClusterAbstract wrappedCluster;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    protected final boolean[] flippedAtoms;
    private IVector childAtomVector;
}
