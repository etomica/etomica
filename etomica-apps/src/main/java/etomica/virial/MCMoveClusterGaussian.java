package etomica.virial;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.potential.PotentialMaster;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.util.random.IRandom;

public class MCMoveClusterGaussian extends MCMoveAtom {
    public MCMoveClusterGaussian(IRandom random, Space space, double [][] meanPosition, double [][] standardDeviation) {
        super(random, null, space);
        this.standardDeviation = standardDeviation;
        this.meanPosition = meanPosition;
    }

    public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        IAtomList leafAtoms = box.getLeafList();
        for(int i = 0; i <leafAtoms.size(); i++) {
            Vector r = leafAtoms.get(i).getPosition();
            for(int j = 0; j < box.getSpace().D(); j++){
                r.setX(j, meanPosition[i][j] + random.nextGaussian() * standardDeviation[i][j]);
            }
            if (imposePBC) r.PE(box.getBoundary().centralImage(r));
        }
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }

    public void acceptNotify() {
        super.acceptNotify();
        ((BoxCluster)box).acceptNotify();
    }

    public double getChi(double temperature) {
        return (uOld==0.0) ? Double.POSITIVE_INFINITY : uNew/uOld;
    }

    public void setDoImposePBC(boolean doImposePBC) {
        imposePBC = doImposePBC;
    }

    public void rejectNotify() {
    }

    protected boolean imposePBC = false;
    protected double [][] standardDeviation;
    protected double [][] meanPosition;
}
