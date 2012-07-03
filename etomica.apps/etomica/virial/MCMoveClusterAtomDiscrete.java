package etomica.virial;

import etomica.api.IRandom;
import etomica.api.IVectorMutable;
import etomica.integrator.mcmove.MCMoveAtom;
import etomica.space.ISpace;

/**
 * Extension of MCMoveClusterAtom that moves only the second atom and only
 * along the x axis.  The atom is moved in discrete steps such that its
 * position is always an integer multiple of the discrete step size.
 * 
 * @author Andrew Schultz
 */
public class MCMoveClusterAtomDiscrete extends MCMoveAtom {

    public MCMoveClusterAtomDiscrete(IRandom random, ISpace _space, double dr) {
        super(random, null, _space);
        setStepSize(1.2);
        setStepSizeMin(2*dr);
        this.dr = dr;
	}

    public boolean doTrial() {
        uOld = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        int imax = (int)Math.ceil(stepSize / dr);
        int idr = random.nextInt(2*imax) - imax;
        if (idr >= 0) idr++;
        IVectorMutable p1 = box.getLeafList().getAtom(1).getPosition();
        oldR = p1.getX(0);
        int iOldR = (int)Math.round(oldR/dr);
        newR = (iOldR + idr)*dr;
        p1.setX(0, newR);
        ((BoxCluster)box).trialNotify();
        uNew = ((BoxCluster)box).getSampleCluster().value((BoxCluster)box);
        return true;
    }
    
    public double getA() {
        if (uOld == 0) return Double.POSITIVE_INFINITY;
        double ratio = uNew/uOld;
        // we need to give r=0 double weight since we visit r=-1 and r=+1 
        if (oldR == 0) ratio *= 0.5;
        else if (newR == 0) ratio *= 2;
        return ratio;
    }

    public double getB() {
        return 0;
    }
    
    public void rejectNotify() {
        IVectorMutable p1 = box.getLeafList().getAtom(1).getPosition();
        p1.setX(0, oldR);
        ((BoxCluster)box).rejectNotify();
    }
    
    public void acceptNotify() {
    	((BoxCluster)box).acceptNotify();
    }

    private static final long serialVersionUID = 1L;

    protected double dr, oldR, newR;
}
