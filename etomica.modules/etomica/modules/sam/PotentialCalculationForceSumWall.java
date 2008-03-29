package etomica.modules.sam;

import etomica.api.IAtomSet;
import etomica.api.IPotential;
import etomica.api.IVector;
import etomica.integrator.IntegratorBox;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialSoft;

public class PotentialCalculationForceSumWall extends
        PotentialCalculationForceSum {

    public PotentialCalculationForceSumWall(P1WCAWall wallPotential) {
        this.wallPotential = wallPotential;
    }
    
    public P1WCAWall getWallPotential() {
        return wallPotential;
    }

    public void reset() {
        super.reset();
        wallForceSum = 0;
    }
    
    public double getWallForce() {
        return wallForceSum;
    }
    
    public void doCalculation(IAtomSet atoms, IPotential potential) {
        PotentialSoft potentialSoft = (PotentialSoft)potential;
        int nBody = potential.nBody();
        IVector[] f = potentialSoft.gradient(atoms);
        switch(nBody) {
            case 1:
                ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(f[0]);
                if (potential == wallPotential) {
                    wallForceSum += f[0].x(wallPotential.getWallDim());
                }
                break;
            case 2:
                ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(0))).force().ME(f[0]);
                ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(1))).force().ME(f[1]);
                break;
            default:
                //XXX atoms.count might not equal f.length.  The potential might size its 
                //array of vectors to be large enough for one IAtomSet and then not resize it
                //back down for another IAtomSet with fewer atoms.
                for (int i=0; i<atoms.getAtomCount(); i++) {
                    ((IntegratorBox.Forcible)integratorAgentManager.getAgent(atoms.getAtom(i))).force().ME(f[i]);
                }
        }
    }

    protected final P1WCAWall wallPotential;
    protected double wallForceSum;
}
