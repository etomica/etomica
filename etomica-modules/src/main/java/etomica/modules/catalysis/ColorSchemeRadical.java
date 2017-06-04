/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomType;
import etomica.simulation.Simulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeAgentManager;
import etomica.graphics.ColorSchemeByType;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;

public class ColorSchemeRadical extends ColorSchemeByType {

    public ColorSchemeRadical(Simulation sim, AtomLeafAgentManager agentManager) {
        super(sim);
        this.agentManager = agentManager;
        radicalColorMap = new AtomTypeAgentManager(this, sim);
        fullBondColorMap = new AtomTypeAgentManager(this, sim);
    }

    public Color getAtomColor(IAtom atom) {
        CatalysisAgent agent = (CatalysisAgent)agentManager.getAgent(atom);
        if (agent != null) {
            if (agent.isRadical) {
                return (Color)radicalColorMap.getAgent(atom.getType());
            }
            else if (agent.bondedAtom2 != null) {
                return (Color)fullBondColorMap.getAgent(atom.getType());
            }
        }
        return super.getAtomColor(atom);
    }

    public void setRadicalColor(IAtomType type, Color color) {
        radicalColorMap.setAgent(type, color);
    }
    
    public void setFullBondColor(IAtomType type, Color color) {
        fullBondColorMap.setAgent(type, color);
    }

    private static final long serialVersionUID = 1L;
    protected final AtomLeafAgentManager agentManager;
    protected final AtomTypeAgentManager radicalColorMap, fullBondColorMap;
}
