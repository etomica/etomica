/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.potential.P2HardGeneric;
import etomica.space.Boundary;
import etomica.space.Space;
import etomica.space.Vector;
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
public class P2SquareWellBondingCO extends P2HardGeneric {

    protected final Boundary boundary;

    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected final Vector dvOO, drOO;
    protected final AtomLeafAgentManager<CatalysisAgent> agentManager;

    protected int nSurfaceSites;
    protected double epsilonBarrier;
    protected double epsilonBonding;
    protected final double minOOr2;
    protected P2SquareWellSurface potentialCS;

    public P2SquareWellBondingCO(Space space, AtomLeafAgentManager<CatalysisAgent> agentManager, double coreDiameter, double lambda, double epsilon,
                                 int nSurfaceSites, double epsilonBarrier, double epsilonBonding, double minOOr, Boundary boundary) {
        super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon});
        this.agentManager = agentManager;
        this.boundary = boundary;
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        this.nSurfaceSites = nSurfaceSites;
        this.epsilonBarrier = epsilonBarrier;
        this.epsilonBonding = epsilonBonding;
        this.minOOr2 = minOOr * minOOr;
        drOO = space.makeVector();
        dvOO = space.makeVector();
    }

    public double getRange() {
        return wellDiameter;
    }

    @Override
    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {
        CatalysisAgent agent1 = agentManager.getAgent(atom1);
        CatalysisAgent agent2 = agentManager.getAgent(atom2);

        int newState = oldState;
        if (oldState == 1 && bij < 0) {
            virial[0] = 2.0 * reducedMass * bij;
        } else {    // Well collision
            if (bij > 0.0) {         // Separating
                double barrier = epsilon;
                double de = epsilon;
                boolean isBondEvent = false;
                if (agent1.nSurfaceBonds > nSurfaceSites &&
                        agent2.nSurfaceBonds > nSurfaceSites) {
                    barrier = epsilonBonding + epsilonBarrier;
                    de = epsilonBonding;
                    isBondEvent = true;
                } else if (agent1.bondedAtom1 == atom2 || agent1.bondedAtom2 == atom2) {
                    barrier = Double.POSITIVE_INFINITY;
                }

                if (ke < barrier) {     // Not enough kinetic energy to escape
                    virial[0] = 2.0 * reducedMass * bij;
                } else {                 // Escape
                    virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * epsilon / reducedMass));
                    du[0] = de;
                    if (isBondEvent) {
                        if (agent1.bondedAtom1 == atom1) {
                            agent1.bondedAtom1 = agent1.bondedAtom2;
                            agent1.bondedAtom2 = null;
                        } else {
                            agent1.bondedAtom2 = null;
                        }
                        agent1.isRadical = true;
                        if (agent2.bondedAtom1 == atom1) {
                            agent1.bondedAtom1 = agent2.bondedAtom2;
                            agent2.bondedAtom2 = null;
                        } else {
                            agent2.bondedAtom2 = null;
                        }
                        agent2.isRadical = true;
                    }
                    newState = oldState + 1;
                }
            } else {
                double barrier = -epsilon;
                double de = -epsilon;
                boolean isBondEvent = false;
                if (agent1.isRadical && agent2.isRadical) {
                    IAtomKinetic unbondedO = atom1;
                    IAtomKinetic alreadyBondedO = (IAtomKinetic) agent2.bondedAtom1;
                    int nCSites = agent2.nSurfaceBonds;
                    if (alreadyBondedO == null) {
                        alreadyBondedO = (IAtomKinetic) agent1.bondedAtom1;
                        nCSites = agent1.nSurfaceBonds;
                        unbondedO = atom2;
                    }
                    dvOO.Ev1Mv2(unbondedO.getVelocity(), alreadyBondedO.getVelocity());

                    drOO.Ev1Mv2(unbondedO.getPosition(), alreadyBondedO.getPosition());
                    drOO.PEa1Tv1(falseTime, dvOO);
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
                if (ke > barrier) {   // Approach/capture
                    virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * de / reducedMass));
                    du[0] = de;
                    if (isBondEvent) {
                        if (agent1.bondedAtom1 == null) {
                            agent1.bondedAtom1 = atom2;
                        } else {
                            agent1.bondedAtom2 = atom2;
                        }
                        if (agent2.bondedAtom1 == null) {
                            agent2.bondedAtom1 = atom1;
                        } else {
                            agent2.bondedAtom2 = atom1;
                        }
                        agent1.isRadical = false;
                        agent2.isRadical = false;
                    }
                    newState = oldState - 1;
                } else {                     // Not enough kinetic energy to overcome square-shoulder
                    virial[0] = 2.0 * reducedMass * bij;
                }
            }
        }
        return newState;
    }

    public double u(double r2) {
        throw new RuntimeException("Need to call the non-spherical method");
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        CatalysisAgent agent0 = agentManager.getAgent(atom1);
        if (agent0.bondedAtom1 == atom2 || agent0.bondedAtom2 == atom2) {
            return -epsilonBonding;
        }
        return super.u(dr12, atom1, atom2);
    }

    /**
     * Accessor method for core diameter.
     */
    public double getCoreDiameter() {
        return coreDiameter;
    }

    /**
     * Accessor method for core diameter.
     * Well diameter is defined as a multiple (lambda) of this, and is updated when core diameter is changed
     */
    public void setCoreDiameter(double c) {
        if (c < 0) {
            throw new IllegalArgumentException("diameter must not be negative");
        }
        coreDiameter = c;
        coreDiameterSquared = c * c;
        wellDiameter = coreDiameter * lambda;
        wellDiameterSquared = wellDiameter * wellDiameter;
    }

    public Dimension getCoreDiameterDimension() {
        return Length.DIMENSION;
    }

    /**
     * Accessor method for well-diameter multiplier.
     */
    public double getLambda() {
        return lambda;
    }

    /**
     * Accessor method for well-diameter multiplier.
     * Well diameter is defined as this multiple of core diameter, and is updated when
     * this is changed
     */
    public void setLambda(double lam) {
        if (lam <= 1.0) throw new IllegalArgumentException("Square-well lambda must be greater than 1.0");
        lambda = lam;
        wellDiameter = coreDiameter * lambda;
        wellDiameterSquared = wellDiameter * wellDiameter;
    }

    public Dimension getLambdaDimension() {
        return Null.DIMENSION;
    }

    /**
     * Accessor method for depth of well
     */
    public double getEpsilon() {
        return epsilon;
    }

    /**
     * Accessor method for depth of well
     */
    public void setEpsilon(double eps) {
        epsilon = eps;
    }

    public Dimension getEpsilonDimension() {
        return Energy.DIMENSION;
    }

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
  
