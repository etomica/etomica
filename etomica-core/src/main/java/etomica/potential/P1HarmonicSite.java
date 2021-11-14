/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomList;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.units.dimensions.CompoundDimension;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;

/**
 * Potential in which attaches a harmonic spring between each affected atom a
 * nominal site for that atom, as defined by an agent manager.
 *
 * @author David Kofke
 * @author Andrew Schultz
 */
 
public class P1HarmonicSite extends Potential1 implements PotentialSoft, IPotentialField {

    private final Space space;
    private double w = 100.0;
    private final Vector[] force;
    protected final AtomLeafAgentManager<? extends Vector> atomAgentManager;
    
    public P1HarmonicSite(Space space, AtomLeafAgentManager<? extends Vector> agentManager) {
        super(space);
        this.space = space;
        force = new Vector[]{space.makeVector()};
        atomAgentManager = agentManager;
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
        Vector x0 = atomAgentManager.getAgent(a.get(0));
        return w*a.get(0).getPosition().Mv1Squared(x0);
    }

    public double u(IAtom a) {
        Vector x0 = atomAgentManager.getAgent(a);
        return w*a.getPosition().Mv1Squared(x0);
    }

    @Override
    public double udu(IAtom a, Vector f) {
        Vector r = a.getPosition();
        Vector x0 = atomAgentManager.getAgent(a);
        Vector dr = space.makeVector();
        dr.Ev1Mv2(r,x0);
        f.PEa1Tv1(-2*w, dr);
        return w*dr.squared();
    }

    private Vector[] gradient(IAtomList a){
        Vector r = a.get(0).getPosition();
        Vector x0 = atomAgentManager.getAgent(a.get(0));
        force[0].Ev1Mv2(r,x0);
        force[0].TE(2*w);
            
        return force;
    }

}
   
