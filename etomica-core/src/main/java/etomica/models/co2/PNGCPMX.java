/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.co2;

import etomica.atom.AtomType;
import etomica.space.Space;

import java.util.Map;

/**
 * GCPM potential that allows combining parameters to be adjusted.
 * All combining rules are adjusted the same way (not what you want for a
 * ternary mixture).
 * 
 * @author Andrew Schultz
 */
public class PNGCPMX extends PNGCPM {

    protected final double kSigma, kEpsilon, kGamma;
    
    public PNGCPMX(Space space, Map<AtomType, GCPMAgent> typeManager,
                   int nAtomTypes, double kSigma, double kEpsilon, double kGamma) {
        super(space, typeManager, nAtomTypes);
        this.kSigma = kSigma;
        this.kEpsilon = kEpsilon;
        this.kGamma = kGamma;
    }

    public GCPMAgent getPairAgent(AtomType type1, AtomType type2) {
        int idx1 = type1.getIndex();
        int idx2 = type2.getIndex();
        if (pairAgents[idx1][idx2] != null) return pairAgents[idx1][idx2];
        GCPMAgent agent1 = typeManager.get(type1);
        if (idx1==idx2) {
            pairAgents[idx1][idx2] = agent1;
            return agent1;
        }
        GCPMAgent agent2 = typeManager.get(type2);
        double sigma = kSigma*0.5*(agent1.sigma + agent2.sigma);
        double epsilon = kEpsilon*2*agent1.epsilon*agent2.epsilon;
        if (epsilon>0) epsilon /= (agent1.epsilon + agent2.epsilon);
        double gamma = kGamma*0.5*(agent1.gamma + agent2.gamma);
        double tau = Math.sqrt(0.5*(agent1.tau*agent1.tau + agent2.tau*agent2.tau));
        pairAgents[idx1][idx2] = new GCPMAgent(sigma, epsilon, tau, gamma, agent1.charge, agent2.charge, 0, 0, 0);
        pairAgents[idx2][idx1] = pairAgents[idx1][idx2];
        return pairAgents[idx1][idx2];
    }
}
