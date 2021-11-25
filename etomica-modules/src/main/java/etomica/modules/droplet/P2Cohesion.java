/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.droplet;

import etomica.atom.AtomTest;
import etomica.atom.IAtom;
import etomica.potential.Potential2Soft;
import etomica.space.Vector;

import java.util.function.Predicate;

/**
 * Cohesive potential for mesoscale droplet simulation
 * @author Andrew Schultz
 */
public class P2Cohesion implements Potential2Soft {

    @Override
    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        if (useSurfaceOnly && (liquidFilter.test(atom1) || liquidFilter.test(atom2))) {
            return 0;
        }
        return u(dr12.squared());
    }

    public double d2u(double r2) {
        if (r2 > epsilonSq) {
            return 0;
        }
        return r2*fac*(1-3*r2/epsilonSq)*dv;
    }

    public double du(double r2) {
        if (r2 > epsilonSq) {
            return 0;
        }
        return r2*fac*(1-r2/epsilonSq)*dv;
    }

    public double u(double r2) {
        if (r2 > epsilonSq) {
            return 0;
        }
        return 0.5*r2*fac*(1-0.5*r2/epsilonSq)*dv;
    }
    
    public double getRange() {
        return epsilon;
    }

    public void setEpsilon(double newEpsilon) {
        if (newEpsilon < 0) {
            throw new RuntimeException("Ooops");
        }
        epsilon = newEpsilon;
        epsilonSq = epsilon*epsilon;
        fac = 192/Math.PI/(epsilonSq*epsilonSq*epsilonSq);
    }
    
    public double getEpsilon() {
        return epsilon;
    }
    
    public void setDv(double newDv) {
        dv = newDv;
    }
    
    public double getDv() {
        return dv;
    }

    public void setLiquidFilter(AtomTest newLiquidFilter) {
        liquidFilter = newLiquidFilter;
    }

    public Predicate<etomica.atom.IAtom> getLiquidFilter() {
        return liquidFilter;
    }

    public void setUseSurfaceOnly(boolean newUseSurfaceOnly) {
        useSurfaceOnly = newUseSurfaceOnly;
    }
    
    public boolean getUseSurfaceOnly() {
        return useSurfaceOnly;
    }
    
    protected double epsilon, epsilonSq;
    protected double fac, dv;
    protected boolean useSurfaceOnly;
    protected AtomTest liquidFilter;
}
