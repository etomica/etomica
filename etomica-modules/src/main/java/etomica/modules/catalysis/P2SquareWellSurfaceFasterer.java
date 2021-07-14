/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.IAtom;
import etomica.atom.IAtomKinetic;
import etomica.atom.IAtomList;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.potential.P2HardGeneric;
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
public class P2SquareWellSurfaceFasterer extends P2HardGeneric {

    protected double coreDiameter, coreDiameterSquared;
    protected double wellDiameter, wellDiameterSquared;
    protected double lambda; //wellDiameter = coreDiameter * lambda
    protected double epsilon;
    protected final AtomLeafAgentManager<CatalysisAgent> agentManager;
    protected int minRadicalSites;

    public P2SquareWellSurfaceFasterer(AtomLeafAgentManager<CatalysisAgent> agentManager, double coreDiameter, double lambda, double epsilon, int minRadicalSites) {
        super(new double[]{coreDiameter, coreDiameter * lambda}, new double[]{Double.POSITIVE_INFINITY, -epsilon});
        this.agentManager = agentManager;
        setCoreDiameter(coreDiameter);
        setLambda(lambda);
        setEpsilon(epsilon);
        this.minRadicalSites = minRadicalSites;
    }

    public double getRange() {
        return wellDiameter;
    }

    @Override
    protected int decideBump(IAtomKinetic atom1, IAtomKinetic atom2, int oldState, boolean core, double ke, double reducedMass, double bij, double r2, double[] du, double[] virial, double falseTime) {

        double rm0 = atom1.getType().rm();
        CatalysisAgent agent;
        if (rm0 == 0) {
            agent = agentManager.getAgent(atom2);
        } else {
            agent = agentManager.getAgent(atom1);
        }
        int newState = oldState;
        if (oldState == 1 && bij < 0) {
            virial[0] = 2.0 * reducedMass * bij;
        } else {    // Well collision
            // ke is kinetic energy due to components of velocity
            if (bij > 0.0) {         // Separating
                if (!agent.isRadical && agent.bondedAtom2 != null) {
                    // turn off C-surface attraction for OCO
                    newState = oldState + 1;
                } else if (ke < epsilon || (agent.isRadical && agent.nSurfaceBonds <= minRadicalSites)) {     // Not enough kinetic energy to escape
                    virial[0] = 2.0 * reducedMass * bij;
                } else {                 // Escape
                    du[0] = epsilon;
                    virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * epsilon / reducedMass));
                    if (Double.isNaN(virial[0])) throw new RuntimeException("NaN virial");
                    newState = oldState + 1;
                    agent.nSurfaceBonds--;
                }
            } else {   // Approach/capture
                if (!agent.isRadical && agent.bondedAtom2 != null) {
                    // turn off C-surface attraction for OCO
                } else {
                    du[0] = -epsilon;
                    virial[0] = reducedMass * (bij + (core ? +1 : -1) * Math.sqrt(bij * bij - 2.0 * r2 * (-epsilon) / reducedMass));
                    if (Double.isNaN(virial[0])) throw new RuntimeException("NaN virial");
                    agent.nSurfaceBonds++;
                }
                newState = oldState - 1;
            }
        }
        return newState;
    }

    public double energy(IAtomList pair) {
        IAtom atom0 = pair.get(0);
        IAtom atom1 = pair.get(1);
        double rm0 = atom0.getType().rm();
        CatalysisAgent agent;
        if (rm0 == 0) {
            agent = agentManager.getAgent(atom1);
        } else {
            agent = agentManager.getAgent(atom0);
        }
        if (!agent.isRadical && agent.bondedAtom2 != null) {
            // turn off C-surface attraction for OCO
            return 0;
        }
        return super.energy(pair);
    }

    public double u(Vector dr12, IAtom atom1, IAtom atom2) {
        double rm0 = atom1.getType().rm();
        CatalysisAgent agent;
        if (rm0 == 0) {
            agent = agentManager.getAgent(atom2);
        } else {
            agent = agentManager.getAgent(atom1);
        }
        if (!agent.isRadical && agent.bondedAtom2 != null) {
            // turn off C-surface attraction for OCO
            return 0;
        }
        return super.u(dr12, atom1, atom2);
    }

    /**
     * Returns infinity if overlapping, -epsilon if otherwise less than well diameter, or zero if neither.
     */
    public double u(double r2) {
        throw new RuntimeException("Need to call the non-spherical method");
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

    public int getMinRadicalSites() {
        return minRadicalSites;
    }

    public void setMinRadicalSites(int newMinRadicalSites) {
        if (newMinRadicalSites < 1) {
            throw new RuntimeException("Must be positive");
        }
        minRadicalSites = newMinRadicalSites;
    }
}
  
