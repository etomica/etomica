/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.potential;

import etomica.api.IAtom;
import etomica.api.IAtomKinetic;
import etomica.api.IAtomList;
import etomica.api.IBox;
import etomica.api.IPotentialAtomic;
import etomica.integrator.IntegratorBox;
import etomica.space.ISpace;
import etomica.space.Tensor;

/**
 * Calculates the pressure tensor by calculating the force on each atom, along
 * with including the kinetic portion (from the velocities or an Integrator).
 * If the simulation is non-dynamic (MC), the Integrator must be provided.
 */
public class PotentialCalculationPressureTensor implements PotentialCalculation {

    protected final Tensor pressureTensor;
    protected final Tensor workTensor;
    protected final ISpace space;
    protected IAtomList leafList;
    protected IntegratorBox integrator;
    protected boolean warningPrinted;
    protected double temperature;
    protected final Tensor I;
    protected boolean doNonEquilibrium;
    
    public PotentialCalculationPressureTensor(ISpace space) {
        this.space = space;
        pressureTensor = space.makeTensor();
        workTensor = space.makeTensor();
        I = space.makeTensor();
        for (int i=0; i<space.D(); i++) {
            I.setComponent(i,i,1);
        }
    }

    /**
     * @param doNonEquilibrium
     *      If true, the kinetic contribution to the pressure tensor will be
     *      computed from atomic velocities rather than the temperature.
     */
    public void setDoNonEquilibrium(boolean doNonEquilibrium) {
        this.doNonEquilibrium = doNonEquilibrium;
    }

    /**
	 * Adds the pressure tensor contribution based on the forces acting on each
     * pair of atoms produced by the iterator.
	 */
	public void doCalculation(IAtomList atoms, IPotentialAtomic potential) {
		((PotentialSoft)potential).gradient(atoms, pressureTensor);
	}

    public void setBox(IBox newBox) {
        leafList = newBox.getLeafList();
    }

    public void zeroSum() {
        pressureTensor.E(0);
    }

    public void setTemperature(double newTemperature) {
        this.temperature = newTemperature;
    }

    /**
     * Sets an integrator to use a source for the temperature to compute the
     * kinetic portion of the pressure.  If running a dynamic simulation
     * (where the Atoms have velocities), this method should not be called.
     */
    public void setIntegrator(IntegratorBox newIntegrator) {
        integrator = newIntegrator;
    }

    /**
     * Returns the pressure tensor based on a previous call to 
     * PotentialMaster.calculate
     */
    public Tensor getPressureTensor() {
        if (leafList.getAtomCount() == 0) {
            return pressureTensor;
        }

        if (doNonEquilibrium) {
            
            // use the velocities
            int nLeaf = leafList.getAtomCount();
            for (int iLeaf=0; iLeaf<nLeaf; iLeaf++) {
                IAtomKinetic atom = (IAtomKinetic)leafList.getAtom(iLeaf);
                workTensor.Ev1v2(atom.getVelocity(), atom.getVelocity());
                workTensor.TE(((IAtom)atom).getType().getMass());
                pressureTensor.PE(workTensor);
            }
            return pressureTensor;
        }

        // or just include ideal gas term
        double T = integrator != null ? integrator.getTemperature() : temperature;
        pressureTensor.PEa1Tt1(leafList.getAtomCount()*T,I);
        return pressureTensor;
    }
}
