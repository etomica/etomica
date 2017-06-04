/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomType;
import etomica.simulation.Simulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeAgentManager;
import etomica.graphics.ColorScheme;
import etomica.util.Arrays;

/**
 * Color scheme for stepwise growth, based on the AtomType (alcohol vs. acid)
 * and the number of bonds the atom can form.
 * 
 * @author Andrew Schultz
 */
public class ColorSchemeStepWise extends ColorScheme implements AtomTypeAgentManager.AgentSource {

    public ColorSchemeStepWise(Simulation sim, AtomLeafAgentManager bondingAgentManager) {
        super();
        simulation = sim;
        colorMaps = new AtomTypeAgentManager[0];
        this.bondingAgentManager = bondingAgentManager;
    }
    
    public Color getAtomColor(IAtom atom) {
        IAtom[] nbrs = (IAtom[])bondingAgentManager.getAgent(atom);
        if (nbrs != null && colorMaps.length > nbrs.length) {
            return (Color)colorMaps[nbrs.length].getAgent(atom.getType());
        }
        // we weren't told how to deal with any atom type with this many bonds.
        return ColorScheme.DEFAULT_ATOM_COLOR;
    }

    /**
     * Sets atoms of the given type and number of bonds to be the given color.
     */
    public void setColor(IAtomType type, int nBonds, Color color) {
        if (nBonds >= colorMaps.length) {
            int oldLength = colorMaps.length;
            colorMaps = (AtomTypeAgentManager[])Arrays.resizeArray(colorMaps, nBonds+1);
            for (int i=oldLength; i<colorMaps.length; i++) {
                colorMaps[i] = new AtomTypeAgentManager(this, simulation);
            }
        }
        colorMaps[nBonds].setAgent(type, color);
    }
    
    public Class getSpeciesAgentClass() {
        return Color.class;
    }

    public Object makeAgent(IAtomType type) {
        return ColorScheme.DEFAULT_ATOM_COLOR;
    }

    public void releaseAgent(Object agent, IAtomType type) {}

    protected AtomTypeAgentManager[] colorMaps;
    protected final AtomLeafAgentManager bondingAgentManager;
    protected final Simulation simulation;
}
