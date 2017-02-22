/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.data.meter;

import etomica.api.IAtom;
import etomica.api.IAtomKinetic;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IMoleculeList;
import etomica.api.ISimulation;
import etomica.api.IVector;
import etomica.api.IVectorMutable;
import etomica.atom.IMoleculeKinetic;
import etomica.atom.IMoleculeOrientedKinetic;
import etomica.data.DataSourceScalar;
import etomica.space.ISpace;
import etomica.space3d.IOrientationFull3D;
import etomica.space3d.RotationTensor3D;
import etomica.species.ISpeciesOriented;
import etomica.units.Energy;

/**
 * Meter for the total kinetic energy in a box
 * Computes total KE by summing values of KE returned by every atom in the box.
 * Kinetic energy for rigid molecules is calculated based on the its
 * translational and rotational energy.
 * 
 * @author Andrew Schultz
 */
public class MeterKineticEnergyRigid extends DataSourceScalar {

    public MeterKineticEnergyRigid(ISpace space, ISimulation sim) {
        this(space, sim, null);
    }
    
    public MeterKineticEnergyRigid(ISpace space, ISimulation sim, IBox box) {
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
            if (moleculeList.getMoleculeCount() == 0) {
                continue;
            }
            IMolecule molecule0 = moleculeList.getMolecule(0);
            if (molecule0 instanceof IMoleculeOrientedKinetic) {
                for (int j=0; j<moleculeList.getMoleculeCount(); j++) {
                    IMoleculeOrientedKinetic moleculeOrientedKinetic = (IMoleculeOrientedKinetic)moleculeList.getMolecule(j);
                    double mass = ((ISpeciesOriented)molecule0.getType()).getMass();
                    if (Double.isInfinite(mass)) {
                        continue;
                    }
                    ke += 0.5*mass*((IMoleculeKinetic)molecule0).getVelocity().squared();
                    IVector moment = ((ISpeciesOriented)molecule0.getType()).getMomentOfInertia();
        
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
                for (int j=0; j<moleculeList.getMoleculeCount(); j++) {
                    IMolecule molecule = moleculeList.getMolecule(j);
                    IAtomList children = molecule.getChildList();
                    for (int iLeaf=0; iLeaf<children.getAtomCount(); iLeaf++) {
                        IAtomKinetic a = (IAtomKinetic)children.getAtom(iLeaf);
                        double mass = ((IAtom)a).getType().getMass();
                        if(mass == Double.POSITIVE_INFINITY) continue;
        //                    System.out.println("force: "+((MyAgent)a.ia).force.toString());
                        IVectorMutable velocity = a.getVelocity();
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
    public void setBox(IBox newBox) {
        box = newBox;
    }
    
    /**
     * @return Returns the box.
     */
    public IBox getBox() {
        return box;
    }

    private static final long serialVersionUID = 1L;
    protected final ISimulation sim;
    protected IBox box;
    protected final IVectorMutable angularVelocity;
    protected final RotationTensor3D rotationTensor;
 }