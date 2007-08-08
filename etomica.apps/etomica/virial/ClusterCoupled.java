package etomica.virial;

import etomica.atom.AtomSet;
import etomica.atom.IAtomGroup;
import etomica.atom.IAtomPositioned;

public class ClusterCoupled implements ClusterAbstract {

    public ClusterCoupled(ClusterAbstract cluster) {
        wrappedCluster = cluster;
    }

    public ClusterAbstract makeCopy() {
        return new ClusterCoupled(wrappedCluster);
    }

    public int pointCount() {
        return wrappedCluster.pointCount();
    }

    public ClusterAbstract getSubCluster() {
        return wrappedCluster;
    }
    
    public double value(BoxCluster box) {
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
        
        double v1 = wrappedCluster.value(box);
        double r = invert(cPairs);
        double v2;
        if (r == -1) {
            return v1;
        }
        v2 = wrappedCluster.value(box);
        foo += v2/Math.pow(r,3)*v1;
        foo2 += v1*v1;
        double ri = 1/r;
        value = (v1 + v2*(ri*ri*ri))*0.5;
        if (v1*v2 < 0) {
//            System.out.println("r "+r+" v "+v1+" v2 "+v2);
//            System.out.print("v "+v1*(r*r)+" v2 "+v2/(ri*ri));
//            System.out.println(" w "+v+" r "+r);
        }
//        v /= (1 + r2*r2);
        invert(cPairs);
        return value;
    }
    
    private double invert(CoordinatePairSet cPairs) {
        double minmax = -1;
        IAtomPositioned minMaxAtom = null;
        for (int i=1; i<atomList.getAtomCount(); i++) {
            IAtomPositioned atom = (IAtomPositioned)atomList.getAtom(i);
            double r = atom.getPosition().squared(); // sqrt
            if (r == 1.0) {
                return 1.0;
            }
            if (r > 1.0) {
                r = 1/r;
            }
//            if (r > 2.0) {
//                continue;
//            }
//            if (r > 1.0) {
//                r = 2-r;
//            }
            if (r > minmax) {
                minMaxAtom = atom;
                minmax = r;
            }
        }
        if (minMaxAtom == null) {
//            System.out.println("skipping");
            return -1.0;
        }
        double r = minMaxAtom.getPosition().squared();  // sqrt
//        System.out.println(minMaxAtom+" "+r);
//        System.out.println("picked "+minMaxAtom+" at "+minMaxAtom.coord.position()+" (r="+r+")");
        minMaxAtom.getPosition().TE(1/r); //(2/r-1);
//        double rN = Math.sqrt(minMaxAtom.coord.position().squared());
//        System.out.println("now at "+minMaxAtom.coord.position()+" (r="+rN+")");
        cPairs.reset();
        return r;
    }
    
    public void setBox(BoxCluster box) {
        atomList = ((IAtomGroup)box.getSpeciesMaster().getAgentList().getAtom(0)).getChildList();
    }

    public void setTemperature(double temperature) {
        wrappedCluster.setTemperature(temperature);
    }
    
    private final ClusterAbstract wrappedCluster;
    private AtomSet atomList;
    protected int cPairID = -1, lastCPairID = -1;
    protected double value, lastValue;
    public double foo, foo2;
}
