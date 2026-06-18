/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.spin.ising;

import etomica.atom.IAtomList;

/**
 * Magnetic spin potential, with an energy defined by
 * <p>
 * U = -J r1 dot r2
 * <p>
 * where J is a coupling parameter, and r1 and r2 are the vectors given by
 * atom.coord.position. It is expected (but not verified here) that these
 * vectors are normalized to unity, and that the simulation integrator's
 * algorithm enforces this constraint.
 *
 * @author David Kofke
 */
public class P2Spin {

    private double coupling;

    public P2Spin() {
        this(1.0);
    }

    public P2Spin(double coupling) {
        setCoupling(coupling);
    }

    /**
     * Returns the energy for the given pair of atoms.
     *
     * @throws ClassCastException if atoms is not an instance of AtomPair
     */
    public double energy(IAtomList atoms) {
        return -coupling
                * atoms.get(0).getPosition().dot(atoms.get(1).getPosition());
    }

    /**
     * Returns 0, becuase potential operates on a lattice and range
     * should not be needed.  The PotentialMasterSite expects all Potentials
     * to have a range and uses the return value to determine whether or not
     * to use site iteration.
     */
    public double getRange() {
        return 0;
    }

    /**
     * @return Returns the coupling.
     */
    public double getCoupling() {
        return coupling;
    }

    /**
     * @param coupling The coupling to set.
     */
    public void setCoupling(double coupling) {
        this.coupling = coupling;
    }

}
