/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.catalysis;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.IAtom;
import etomica.graphics.ColorSchemeByType;
import etomica.modules.catalysis.InteractionTracker.CatalysisAgent;
import etomica.simulation.Simulation;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

public class ColorSchemeRadical extends ColorSchemeByType {

    protected final AtomLeafAgentManager<CatalysisAgent> agentManager;
    private final Map<AtomType, Color> radicalColorMap, fullBondColorMap;
    
    public ColorSchemeRadical(Simulation sim, AtomLeafAgentManager<CatalysisAgent> agentManager) {
        super(sim);
        this.agentManager = agentManager;
        radicalColorMap = new HashMap<>();
        fullBondColorMap = new HashMap<>();
    }

    public Color getAtomColor(IAtom atom) {
        CatalysisAgent agent = agentManager.getAgent(atom);
        if (agent != null) {
            if (agent.isRadical) {
                return radicalColorMap.get(atom.getType());
            }
            else if (agent.bondedAtom2 != null) {
                return fullBondColorMap.get(atom.getType());
            }
        }
        return super.getAtomColor(atom);
    }

    public void setRadicalColor(AtomType type, Color color) {
        radicalColorMap.put(type, color);
    }

    public void setFullBondColor(AtomType type, Color color) {
        fullBondColorMap.put(type, color);
    }
}
