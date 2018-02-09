/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.sam;

import etomica.atom.IAtomList;
import etomica.potential.IPotentialAtomic;
import etomica.potential.PotentialCalculationForceSum;
import etomica.potential.PotentialSoft;
import etomica.space.Vector;

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
    
    public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
        PotentialSoft potentialSoft = (PotentialSoft)potential;
        int nBody = potential.nBody();
        Vector[] f = potentialSoft.gradient(atoms);
        switch(nBody) {
            case 1:
                if (f[0].squared() > 2.5e9) {
                    double scale = 50000/Math.sqrt(f[0].squared());
                    integratorAgentManager.getAgent(atoms.getAtom(0)).PEa1Tv1(-scale, f[0]);
                    if (potential == wallPotential) {
                        wallForceSum += scale*f[0].getX(wallPotential.getWallDim());
                    }
                }
                else {
                    integratorAgentManager.getAgent(atoms.getAtom(0)).ME(f[0]);
                    if (potential == wallPotential) {
                        wallForceSum += f[0].getX(wallPotential.getWallDim());
                    }
                }
                break;
            case 2:
                if (f[0].squared() > 2.5e9) {
                    double scale = 50000/Math.sqrt(f[0].squared());
                    integratorAgentManager.getAgent(atoms.getAtom(0)).PEa1Tv1(-scale, f[0]);
                    integratorAgentManager.getAgent(atoms.getAtom(1)).PEa1Tv1(-scale, f[1]);
                }
                else {
                    integratorAgentManager.getAgent(atoms.getAtom(0)).ME(f[0]);
                    integratorAgentManager.getAgent(atoms.getAtom(1)).ME(f[1]);
                }
                break;
            default:
                // we hit this for bonding potentials
                for (int i=0; i<atoms.getAtomCount(); i++) {
                    integratorAgentManager.getAgent(atoms.getAtom(i)).ME(f[i]);
                }
        }
    }

    protected final P1WCAWall wallPotential;
    protected double wallForceSum;
}
