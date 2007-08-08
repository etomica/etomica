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
    
    public double value(CoordinatePairSet cPairs, AtomPairSet aPairs) {
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
        
        // Calculate original value (neither molecule flipped)
        double vsum = wrappedCluster.value(cPairs,aPairs);
        value = vsum;
//        if (Math.random() > 0.99) System.out.println("value = " + value);
//        if (true) return value;
//        System.out.println("Original value = " + vsum);
        
        // Flip molecule 0 and calculate value
 //       flip(0);
 //       cPairs.reset();
 //       vsum += wrappedCluster.value(cPairs,aPairs);

        // Flip molecule 1 and calculate value
        flip(1);
        cPairs.reset();
        vsum += wrappedCluster.value(cPairs,aPairs);
        
        // Unflip molecule 0 and calculate value
  //      flip(0);
  //      cPairs.reset();
  //      vsum += wrappedCluster.value(cPairs,aPairs);
        
//        if (cPairs.getCPair(0,1).r2() > 25 && Math.random() > 0.999) System.out.println("Squrared O-O distance = " + cPairs.getCPair(0,1).r2() + ", v1 = " + v1 + ", v2 = " + v2);
        value = vsum/2; //4;
        //if (v1*v2 < 0) {
//            System.out.println("r "+r+" v "+v1+" v2 "+v2);
//            System.out.print("v "+v1*(r*r)+" v2 "+v2/(ri*ri));
//            System.out.println(" w "+v+" r "+r);
        //}
//        v /= (1 + r2*r2);
        
        // Unflip molecule 1 - back to original configuration
        flip(1);
        cPairs.reset();
//        System.out.println("Original value back again? = " + wrappedCluster.value(cPairs,aPairs));
        cPairID = cPairs.getID();
        return value;
    }
    
    private void flip(int moleculeNumber) {
    		IAtomGroup flippedMolecule = (IAtomGroup)atomList.getAtom(moleculeNumber);
    		IVector COM = flippedMolecule.getType().getPositionDefinition().position(flippedMolecule);
    		AtomSet childAtoms = flippedMolecule.getChildList();
    		for (int i = 0; i < childAtoms.getAtomCount(); i++) {
    			childAtomVector.Ea1Tv1(2,COM);
    			childAtomVector.ME(((IAtomPositioned)childAtoms.getAtom(i)).getPosition());
    			((IAtomPositioned)childAtoms.getAtom(i)).getPosition().E(childAtomVector);
    		}
     }
    
    public void setPhase(BoxCluster box) {
        atomList = ((AtomGroup)box.getSpeciesMaster().getAgentList().getAtom(0)).getChildList();
    }

    public void setTemperature(double temperature) {
        wrappedCluster.setTemperature(temperature);
    }
    
    private final ClusterAbstract wrappedCluster;
    private AtomSet atomList;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    public double foo, foo2;
    private IVector childAtomVector;
}
