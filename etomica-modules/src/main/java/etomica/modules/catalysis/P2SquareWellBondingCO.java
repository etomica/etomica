/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.potential.Potential2HardSpherical;
import etomica.space.Vector;
import etomica.space.Space;
import etomica.space.Tensor;
import etomica.units.dimensions.Dimension;
import etomica.units.dimensions.Energy;
import etomica.units.dimensions.Length;
import etomica.units.dimensions.Null;

/**
 * Basic square-well potential.
 * Energy is infinite if spheres overlap, is -epsilon if less than lambda*sigma and not overlapping,
 * and is zero otherwise.  Core diameter describes size of hard core; lambda is multiplier to get range of well.
 * Suitable for use in space of any dimension.
 * Can be used with negative value for epsilon to produce square-shoulder potential. 
 */
public class P2SquareWellBondingCO extends Potential2HardSpherical {

    private static final long serialVersionUID = 1L;
    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected double lastCollisionVirial, lastCollisionVirialr2;
    protected final Tensor lastCollisionVirialTensor;
    protected double lastEnergyChange;
    protected final Vector dv;
    protected final Vector dvOO, drOO;
    protected final AtomLeafAgentManager agentManager;

    protected int nSurfaceSites;
    protected double epsilonBarrier;
    protected double epsilonBonding;
    protected final double minOOr2;
    protected P2SquareWellSurface potentialCS;

    public P2SquareWellBondingCO(Space space, AtomLeafAgentManager agentManager, double coreDiameter, double lambda, double epsilon,
                                 int nSurfaceSites, double epsilonBarrier, double epsilonBonding, double minOOr) {
        super(space);
        this.agentManager = agentManager;
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        this.nSurfaceSites = nSurfaceSites;
        this.epsilonBarrier = epsilonBarrier;
        this.epsilonBonding = epsilonBonding;
        this.minOOr2 = minOOr * minOOr;
        dv = space.makeVector();
        lastCollisionVirialTensor = space.makeTensor();
        drOO = space.makeVector();
        dvOO = space.makeVector();
    }

    public double getRange() {
        return wellDiameter;
    }
    
    /**
     * Implements collision dynamics between two square-well atoms.
     * Includes all possibilities involving collision of hard cores, and collision of wells
     * both approaching and diverging
     */
    public void bump(IAtomList pair, double falseTime) {
        IAtomKinetic atom0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic atom1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(atom1.getVelocity(), atom0.getVelocity());
        
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double eps = 1.0e-10;
        double rm0 = atom0.getType().rm();
        double rm1 = atom1.getType().rm();
        double reduced_m = 1.0/(rm0+rm1);
        double nudge = 0;
        CatalysisAgent agent0 = (CatalysisAgent)agentManager.getAgent(atom0);
        CatalysisAgent agent1 = (CatalysisAgent)agentManager.getAgent(atom1);

        if(2*r2 < (coreDiameterSquared+wellDiameterSquared)) {   // Hard-core collision
            lastCollisionVirial = 2.0*reduced_m*bij;
            lastEnergyChange = 0.0;
        }
        else {    // Well collision
            // ke is kinetic energy due to components of velocity
            double ke = bij*bij*reduced_m/(2.0*r2);
            if(bij > 0.0) {         // Separating
                double barrier = epsilon;
                double de = epsilon;
                boolean isBondEvent = false;
                if (agent0.nSurfaceBonds > nSurfaceSites &&
                    agent1.nSurfaceBonds > nSurfaceSites) {
                    barrier = epsilonBonding + epsilonBarrier;
                    de = epsilonBonding;
                    isBondEvent = true;
                }
                else if (agent0.bondedAtom1 == atom1 || agent0.bondedAtom2 == atom1) {
                    barrier = Double.POSITIVE_INFINITY;
                }

                if (ke < barrier) {     // Not enough kinetic energy to escape
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = -eps;
                    lastEnergyChange = 0.0;
                }
                else {                 // Escape
                    lastCollisionVirial = reduced_m*(bij - Math.sqrt(bij*bij - 2.0*r2*de/reduced_m));
                    nudge = eps;
                    lastEnergyChange = de;
                    if (isBondEvent) {
                        if (agent0.bondedAtom1 == atom1) {
                            agent0.bondedAtom1 = agent0.bondedAtom2;
                            agent0.bondedAtom2 = null;
                        }
                        else {
                            agent0.bondedAtom2 = null;
                        }
                        agent0.isRadical = true;
                        if (agent1.bondedAtom1 == atom0) {
                            agent1.bondedAtom1 = agent1.bondedAtom2;
                            agent1.bondedAtom2 = null;
                        }
                        else {
                            agent1.bondedAtom2 = null;
                        }
                        agent1.isRadical = true;
                    }
                }
            }
            else {
                double barrier = -epsilon;
                double de = -epsilon;
                boolean isBondEvent = false;
                if (agent0.isRadical && agent1.isRadical) {
                    IAtomKinetic unbondedO = atom0;
                    IAtomKinetic alreadyBondedO = (IAtomKinetic)agent1.bondedAtom1;
                    int nCSites = agent1.nSurfaceBonds;
                    if (alreadyBondedO == null) {
                        alreadyBondedO = (IAtomKinetic)agent0.bondedAtom1;
                        nCSites = agent0.nSurfaceBonds;
                        unbondedO = atom1;
                    }
                    dvOO.Ev1Mv2(unbondedO.getVelocity(), alreadyBondedO.getVelocity());
                    
                    drOO.Ev1Mv2(unbondedO.getPosition(), alreadyBondedO.getPosition());
                    drOO.PEa1Tv1(falseTime,dvOO);
                    boundary.nearestImage(drOO);

                    double OOr2 = drOO.squared();
                    if (OOr2 > minOOr2) {
                        barrier = epsilonBarrier;
                        // when we form the C-O bond, we turn off the C reactions with the surface
                        // (otherwise the CO2 is too strongly adsorbed)
                        // we need to account for this energy change in the bonding energy here.
                        de = -epsilonBonding + nCSites * potentialCS.getEpsilon();
                        if (de > 0) {
                            barrier += de;
                        }
                        isBondEvent = true;
                    }
                }
                if(ke > barrier) {   // Approach/capture
                    lastCollisionVirial = reduced_m*(bij +Math.sqrt(bij*bij-2.0*r2*de/reduced_m));
                    nudge = -eps;
                    lastEnergyChange = de;
                    if (isBondEvent) {
                        if (agent0.bondedAtom1 == null) {
                            agent0.bondedAtom1 = atom1;
                        }
                        else {
                            agent0.bondedAtom2 = atom1;
                        }
                        if (agent1.bondedAtom1 == null) {
                            agent1.bondedAtom1 = atom0;
                        }
                        else {
                            agent1.bondedAtom2 = atom0;
                        }
                        agent0.isRadical = false;
                        agent1.isRadical = false;
                    }
                }
                else {                     // Not enough kinetic energy to overcome square-shoulder
                    lastCollisionVirial = 2.0*reduced_m*bij;
                    nudge = eps;
                    lastEnergyChange = 0.0;
                }
            }
        }
        lastCollisionVirialr2 = lastCollisionVirial/r2;
        dv.Ea1Tv1(lastCollisionVirialr2,dr);
        atom0.getVelocity().PEa1Tv1( rm0,dv);
        atom1.getVelocity().PEa1Tv1(-rm1,dv);
        atom0.getPosition().PEa1Tv1(-falseTime*rm0,dv);
        atom1.getPosition().PEa1Tv1( falseTime*rm1,dv);
        if(nudge != 0) {
            if (rm0 > 0) {
                atom0.getPosition().PEa1Tv1(-nudge,dr);
            }
            if (rm1 > 0) {
                atom1.getPosition().PEa1Tv1(nudge,dr);
            }
        }
    }//end of bump method

    public double lastCollisionVirial() {
        return lastCollisionVirial;
    }

    public Tensor lastCollisionVirialTensor() {
        lastCollisionVirialTensor.Ev1v2(dr, dr);
        lastCollisionVirialTensor.TE(lastCollisionVirialr2);
        return lastCollisionVirialTensor;
    }

    /**
     * Computes next time of collision of two square-well atoms, assuming free-flight kinematics.
     * Collision may occur when cores collides, or when wells first encounter each other on
     * approach, or when they edge of the wells are reached as atoms diverge.
     */
    public double collisionTime(IAtomList pair, double falseTime) {
        IAtomKinetic coord0 = (IAtomKinetic)pair.get(0);
        IAtomKinetic coord1 = (IAtomKinetic)pair.get(1);
        dv.Ev1Mv2(coord1.getVelocity(), coord0.getVelocity());
        
        dr.Ev1Mv2(coord1.getPosition(), coord0.getPosition());
        dr.PEa1Tv1(falseTime,dv);
        boundary.nearestImage(dr);

        double r2 = dr.squared();
        double bij = dr.dot(dv);
        double v2 = dv.squared();
        double time = Double.POSITIVE_INFINITY;

        if(r2 < wellDiameterSquared) {  // Already inside wells

            if(bij < 0.0) {    // Check for hard-core collision
                double discr = bij*bij - v2 * ( r2 - coreDiameterSquared );
                if(discr > 0) {  // Hard cores collide next
                    time = (-bij - Math.sqrt(discr))/v2;
                }
                else {           // Moving toward each other, but wells collide next
                    discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                    time = (-bij + Math.sqrt(discr))/v2;
                }
            }
            else {           // Moving away from each other, wells collide next
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );  // This is always > 0
                time = (-bij + Math.sqrt(discr))/v2;
            }
        }
        else {              // Outside wells; look for collision at well
            if(bij < 0.0) {
                double discr = bij*bij - v2 * ( r2 - wellDiameterSquared );
                if(discr > 0) {
                    time = (-bij - Math.sqrt(discr))/v2;
                }
            }
        }
        return time + falseTime;
    }

  /**
   * Returns infinity if overlapping, -epsilon if otherwise less than well diameter, or zero if neither.
   */
    public double u(double r2) {
        if (r2 > wellDiameterSquared) return 0.0;
        if (r2 > coreDiameterSquared) return -epsilon;
        return Double.POSITIVE_INFINITY;
    }

    public double energyChange() {return lastEnergyChange;}

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {return coreDiameter;}
    /**
     * Accessor method for core diameter.
     * Well diameter is defined as a multiple (lambda) of this, and is updated when core diameter is changed
     */
    public void setCoreDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        coreDiameter = c;
        coreDiameterSquared = c*c;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getCoreDiameterDimension() {return Length.DIMENSION;}

    /**
     * Accessor method for well-diameter multiplier.
     */
    public double getLambda() {return lambda;}
    /**
     * Accessor method for well-diameter multiplier.
     * Well diameter is defined as this multiple of core diameter, and is updated when 
     * this is changed
     */
    public void setLambda(double lam) {
        if (lam <= 1.0) throw new IllegalArgumentException("Square-well lambda must be greater than 1.0");
        lambda = lam;
        wellDiameter = coreDiameter*lambda;
        wellDiameterSquared = wellDiameter*wellDiameter;
    }
    public Dimension getLambdaDimension() {return Null.DIMENSION;}

    /**
     * Accessor method for depth of well
     */
    public double getEpsilon() {return epsilon;}
    /**
     * Accessor method for depth of well
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }
    public Dimension getEpsilonDimension() {return Energy.DIMENSION;}
    
    public double getBarrier() {
        return epsilonBarrier;
    }

    public void setBarrier(double newBarrier) {
        if (newBarrier < 0) {
            throw new RuntimeException("Must not be negative");
        }
        epsilonBarrier = newBarrier;
    }
    
    public double getEpsilonBonding() {
        return epsilonBonding;
    }

    public void setEpsilonBonding(double newEpsilonBonding) {
        if (newEpsilonBonding < 0) {
            throw new RuntimeException("Must not be negative");
        }
        epsilonBonding = newEpsilonBonding;
    }

    public int getNumSurfaceSites() {
        return nSurfaceSites;
    }

    public void setNumSurfaceSites(int newNumSurfaceSites) {
        if (newNumSurfaceSites < 1) {
            throw new RuntimeException("Must be positive");
        }
        nSurfaceSites = newNumSurfaceSites;
    }

    public void setCSPotential(P2SquareWellSurface newPotentialCS) {
        potentialCS = newPotentialCS;
    }
}
  
