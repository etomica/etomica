package etomica.modules.droplet;

import etomica.api.IPotentialAtomic;
import etomica.potential.Potential2SoftSpherical;
import etomica.space.ISpace;

/**
 * Cohesive potential for mesoscale droplet simulation
 * @author Andrew Schultz
 */
public class P2Cohesion extends Potential2SoftSpherical implements
        IPotentialAtomic {

    public P2Cohesion(ISpace space) {
        super(space);
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

    public double uInt(double rc) {
        return 0;
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
    
    protected double epsilon, epsilonSq;
    protected double fac, dv;
}
