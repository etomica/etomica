/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.box.Box;
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
public class P2SquareWellBondingFasterer extends P2HardGeneric {

    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected final AtomLeafAgentManager<CatalysisAgent> agentManager;
    protected Boundary boundary;

    protected int nSurfaceSites;
    protected double epsilonBarrier;
    protected double epsilonBonding;
    protected final double minOCOr2;

    public P2SquareWellBondingFasterer(Space space, AtomLeafAgentManager<CatalysisAgent> agentManager, double coreDiameter, double lambda, double epsilon,
                                       int nSurfaceSites, double epsilonBarrier, double epsilonBonding, double minOCOr) {
        super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon});
        this.agentManager = agentManager;
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        this.nSurfaceSites = nSurfaceSites;
        this.epsilonBarrier = epsilonBarrier;
        this.epsilonBonding = epsilonBonding;
        minOCOr2 = minOCOr * minOCOr;
    }

    public double getRange() {
        return Math.sqrt(minOCOr2);
    }

    @Override
    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {
        CatalysisAgent agent1 = agentManager.getAgent(atom1);
        CatalysisAgent agent2 = agentManager.getAgent(atom2);
        double barrier, de;
        boolean isBondEvent = false;
        if (oldState == 1 && bij < 0) {
            virial[0] = 2.0 * reducedMass * bij;
            return oldState;
        } else if (bij < 0) {
            de = -epsilon;
            barrier = 0;
            if (agent1.isRadical && agent2.isRadical) {
                barrier = epsilonBarrier;
                de = -epsilonBonding;
                isBondEvent = true;
            } else if (agent1.bondedAtom1 == agent2.bondedAtom1) {
                // OCO, prevent approach
                barrier = Double.POSITIVE_INFINITY;
            }
            if (ke > barrier) {   // Approach/capture
                du[0] = de;
                virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * de / reducedMass));

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
                if (oldState == 1) {
                    throw new RuntimeException("oops");
                }
                return oldState - 1;
            } else {
                virial[0] = 2.0 * reducedMass * bij;
                return oldState;
            }
        } else {
            barrier = -epsilon;
            isBondEvent = false;
            de = -epsilon;
            if (agent1.bondedAtom1 == atom2 || agent1.bondedAtom2 == atom2) {
                // O2
                isBondEvent = true;
                if (agent1.nSurfaceBonds > nSurfaceSites &&
                        agent2.nSurfaceBonds > nSurfaceSites) {
                    barrier = epsilonBonding + epsilonBarrier;
                    de = epsilonBonding;
                } else {
                    // not on the surface, prevent unbonding
                    barrier = Double.POSITIVE_INFINITY;
                }
            }
            if (ke > barrier) {
                if (isBondEvent) {
                    if (agent1.bondedAtom1 == atom2) {
                        agent1.bondedAtom1 = agent1.bondedAtom2;
                        agent1.bondedAtom2 = null;
                    } else {
                        agent1.bondedAtom2 = null;
                    }
                    agent1.isRadical = true;
                    if (agent2.bondedAtom1 == atom1) {
                        agent2.bondedAtom1 = agent2.bondedAtom2;
                        agent2.bondedAtom2 = null;
                    } else {
                        agent2.bondedAtom2 = null;
                    }
                    agent2.isRadical = true;
                }
                du[0] = de;
                virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * de / reducedMass));
                return oldState + 1;
            } else {
                virial[0] = 2.0 * reducedMass * bij;
                return oldState;
            }
        }
    }

    /**
     * Computes next time of collision of two square-well atoms, assuming free-flight kinematics.
     * Collision may occur when cores collides, or when wells first encounter each other on
     * approach, or when they edge of the wells are reached as atoms diverge.
     */
    public double collisionTime(IAtomKinetic atom1, IAtomKinetic atom2, Vector r12, Vector v12, int collisionState) {
        double bij = r12.dot(v12);

        if (bij < 0.0 && collisionState == -1) {
            CatalysisAgent agent0 = agentManager.getAgent(atom1);
            CatalysisAgent agent1 = agentManager.getAgent(atom2);
            if (!agent0.isRadical && !agent1.isRadical && agent0.bondedAtom1 == agent1.bondedAtom1) {
                //OCO
                // Outside wells; look for collision at well
                double r2 = r12.squared();
                double v2 = v12.squared();
                double discr = bij * bij - v2 * (r2 - minOCOr2);
                double time = Double.POSITIVE_INFINITY;
                if (discr > 0) {
                    time = (-bij - Math.sqrt(discr)) / v2;
                }
                return time;
            }
        }
        return super.collisionTime(atom1, atom2, r12, v12, collisionState);
    }

    public double energy(IAtomList pair) {
        IAtom atom0 = pair.get(0);
        IAtom atom1 = pair.get(1);

        Vector dr = Vector.d(atom1.getPosition().getD());
        dr.Ev1Mv2(atom1.getPosition(), atom0.getPosition());
        boundary.nearestImage(dr);
        double r2 = dr.squared();
        if (r2 > wellDiameterSquared) return 0.0;
        if (r2 < coreDiameterSquared) return Double.POSITIVE_INFINITY;
        CatalysisAgent agent0 = agentManager.getAgent(atom0);
        return agent0.bondedAtom1 == atom1 || agent0.bondedAtom2 == atom1 ? -epsilonBonding : -epsilon;
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        int s = getState((IAtomKinetic) atom1, (IAtomKinetic) atom2, dr12);
        if (s != 1) return getEnergyForState(s);
        CatalysisAgent agent0 = agentManager.getAgent(atom1);
        return agent0.bondedAtom1 == atom1 || agent0.bondedAtom2 == atom1 ? -epsilonBonding : -epsilon;
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

    public int getNumSurfaceSites() {
        return nSurfaceSites;
    }

    public void setNumSurfaceSites(int newNumSurfaceSites) {
        if (newNumSurfaceSites < 1) {
            throw new RuntimeException("Must be positive");
        }
        nSurfaceSites = newNumSurfaceSites;
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

    public void setBox(Box box) {
        boundary = box.getBoundary();
    }
}
  
