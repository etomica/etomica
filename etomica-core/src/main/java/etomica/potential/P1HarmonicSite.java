/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.space.Vector;
import etomica.atom.AtomLeafAgentManager;
import etomica.box.BoxAgentManager;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.CompoundDimension;
import etomica.units.Dimension;
import etomica.units.Energy;
import etomica.units.Length;

/**
 * Potential in which attaches a harmonic spring between each affected atom a
 * nominal site for that atom, as defined by an agent manager.
 *
 * @author David Kofke
 * @author Andrew Schultz
 */
 
public class P1HarmonicSite extends Potential1 implements PotentialSoft {
    
    private double w = 100.0;
    private final Vector[] force;
    protected final BoxAgentManager<AtomLeafAgentManager<? extends Vector>> boxAgentManager;
    protected AtomLeafAgentManager<? extends Vector> atomAgentManager;
    
    public P1HarmonicSite(Space space) {
        super(space);
        force = new Vector[]{space.makeVector()};
        boxAgentManager = new BoxAgentManager<AtomLeafAgentManager<? extends Vector>>(null, AtomLeafAgentManager.class);
    }

    public void setAtomAgentManager(Box box, AtomLeafAgentManager<? extends Vector> agentManager) {
        boxAgentManager.setAgent(box, agentManager);
    }
    
    public void setBox(Box box) {
        atomAgentManager = boxAgentManager.getAgent(box);
    }

    public void setSpringConstant(double springConstant) {
        w = springConstant;
    }
    
    public double getSpringConstant() {
        return w;
    }
    
    public Dimension getSpringConstantDimension() {
        return new CompoundDimension(new Dimension[]{Energy.DIMENSION,Length.DIMENSION},new double[]{1,-2});
    }

    public double energy(IAtomList a) {
        Vector x0 = atomAgentManager.getAgent(a.getAtom(0));
        return w*a.getAtom(0).getPosition().Mv1Squared(x0);
    }
    
    public double virial(IAtomList a) {
        return 0.0;
    }

    public Vector[] gradient(IAtomList a){
        Vector r = a.getAtom(0).getPosition();
        Vector x0 = atomAgentManager.getAgent(a.getAtom(0));
        force[0].Ev1Mv2(r,x0);
        force[0].TE(2*w);
            
        return force;
    }
        
    public Vector[] gradient(IAtomList a, Tensor pressureTensor){
        return gradient(a);
    }
}
   
