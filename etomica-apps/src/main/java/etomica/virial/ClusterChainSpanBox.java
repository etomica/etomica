package etomica.virial;

import etomica.atom.IAtomList;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeList;
import etomica.space.Vector;

public class ClusterChainSpanBox implements ClusterAbstract{
    protected final int n;
    protected double[][] pairSigma;
    protected final ClusterAbstract subCluster;

    public ClusterChainSpanBox(int nPoints, double[][] pairSigma, ClusterAbstract subCluster ) {
        this.n = nPoints;
        this.pairSigma = pairSigma;
        this.subCluster = subCluster;
    }

    public ClusterChainSpanBox(int nPoints, double[][] pairSigma) {
        this(nPoints, pairSigma, null);
    }

    @Override
    public ClusterAbstract makeCopy() {
        ClusterChainSpanBox c = new ClusterChainSpanBox(n, pairSigma, subCluster.makeCopy());
        c.setTemperature(1);
        return c;
    }

    @Override
    public int pointCount() {
        return n;
    }

    @Override
    public double value(BoxCluster box) {
        Vector dr = box.space != null ? box.space.makeVector() : null;
        Vector resultant = box.space != null ? box.space.makeVector() : null;
        IAtomList list = box.getLeafList();
        for (int i=0; i<n; i++)
        {
            int j = i == n-1? 0: i+1;
            dr.Ev1Mv2(list.get(j).getPosition(), list.get(i).getPosition());
            double diameter = pairSigma[i][j];
            box.getBoundary().nearestImage(dr);
            if (dr.squared() > diameter * diameter)
                return 0;
            resultant.PE(dr);
        }
//        System.out.println("Resultant"+resultant);
        if (resultant.getX(0) > (box.getBoundary().getBoxSize().getX(0))/2)
            if (subCluster != null) {
                double v = subCluster.value(box);
//                if (v > 0) {
//                    for (int i = 0; i < box.getLeafList().size(); i++) {
//                        System.out.println(box.getLeafList().get(i).getPosition());
//                    }
//                    System.out.println();
//                }
                if (v == 0) {
                    return 1E-200;
                }
                return v;
            }
            else
                return 1;
        else
            return 0;
    }

    @Override
    public void setTemperature(double temperature) {

    }
}
