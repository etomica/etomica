/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
import etomica.data.DataSourceScalar;
import etomica.molecule.IMolecule;
import etomica.molecule.IMoleculeKinetic;
import etomica.molecule.IMoleculeList;
import etomica.molecule.IMoleculeOrientedKinetic;
import etomica.simulation.Simulation;
import etomica.space.Space;
import etomica.space.Vector;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpeciesOriented;
import etomica.units.dimensions.Energy;

/**
 * Meter for the total kinetic energy in a box
 * Computes total KE by summing values of KE returned by every atom in the box.
 * Kinetic energy for rigid molecules is calculated based on the its
 * translational and rotational energy.
 * 
 * @author Andrew Schultz
 */
public class MeterKineticEnergyRigid extends DataSourceScalar {

    public MeterKineticEnergyRigid(Space space, Simulation sim) {
        this(space, sim, null);
    }
    
    public MeterKineticEnergyRigid(Space space, Simulation sim, Box box) {
        super("Kinetic Energy",Energy.DIMENSION);
        angularVelocity = space.makeVector();
        rotationTensor = (RotationTensor3D)space.makeRotationTensor();
        this.box = box;
        this.sim = sim;
    }

	/**
	 * Returns the total kinetic energy summed over all atoms produced by
	 * the iterator when applied to the given box.  Does not include contributions
     * from atoms having infinite mass (it assumes they are stationary).
	 */
    public double getDataAsScalar() {
        if (box == null) throw new IllegalStateException("must call setBox before using meter");
        double ke = 0.0;
        for (int i=0; i<sim.getSpeciesCount(); i++) {
            IMoleculeList moleculeList = box.getMoleculeList(sim.getSpecies(i));
            if (moleculeList.size() == 0) {
                continue;
            }
            IMolecule molecule0 = moleculeList.get(0);
            if (molecule0 instanceof IMoleculeOrientedKinetic) {
                for (int j = 0; j<moleculeList.size(); j++) {
                    IMoleculeOrientedKinetic moleculeOrientedKinetic = (IMoleculeOrientedKinetic)moleculeList.get(j);
                    double mass = ((ISpeciesOriented)molecule0.getType()).getMass();
                    if (Double.isInfinite(mass)) {
                        continue;
                    }
                    ke += 0.5*mass*((IMoleculeKinetic)molecule0).getVelocity().squared();
                    Vector moment = ((ISpeciesOriented)molecule0.getType()).getMomentOfInertia();
        
                    angularVelocity.E(moleculeOrientedKinetic.getAngularVelocity());
                    rotationTensor.setOrientation((IOrientationFull3D)moleculeOrientedKinetic.getOrientation());
                    rotationTensor.transform(angularVelocity);
                    angularVelocity.TE(moment);
        
                    angularVelocity.TE(angularVelocity);
                    angularVelocity.DE(moment);
                    ke += angularVelocity.getX(0) + angularVelocity.getX(1)+ angularVelocity.getX(2);
                }
            }
            else {
                for (int j = 0; j<moleculeList.size(); j++) {
                    IMolecule molecule = moleculeList.get(j);
                    IAtomList children = molecule.getChildList();
                    for (int iLeaf = 0; iLeaf<children.size(); iLeaf++) {
                        IAtomKinetic a = (IAtomKinetic)children.get(iLeaf);
                        double mass = a.getType().getMass();
                        if(mass == Double.POSITIVE_INFINITY) continue;
        //                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                        Vector velocity = a.getVelocity();
                        ke += velocity.squared() * mass;
                    }
                }
            }
        }
        ke *= 0.5;
        return ke;
    }

    /**
     * Sets the box to the given box.
     */
    public void setBox(Box newBox) {
        box = newBox;
    }
    
    /**
     * @return Returns the box.
     */
    public Box getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    protected final Simulation sim;
    protected Box box;
    protected final Vector angularVelocity;
    protected final RotationTensor3D rotationTensor;
 }
