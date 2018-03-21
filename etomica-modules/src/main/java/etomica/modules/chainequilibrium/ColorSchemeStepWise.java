/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomType;
import etomica.atom.IAtom;
import etomica.graphics.ColorScheme;
import etomica.simulation.Simulation;
import etomica.species.ISpecies;

import java.awt.*;

/**
 * Color scheme for stepwise growth, based on the AtomType (alcohol vs. acid) and the number of bonds the atom can
 * form.
 *
 * @author Andrew Schultz
 */
public class ColorSchemeStepWise extends ColorScheme {

    protected final AtomLeafAgentManager<IAtom[]> bondingAgentManager;
    protected final Simulation simulation;
    protected Color[][] colorMaps;

    public ColorSchemeStepWise(Simulation sim, AtomLeafAgentManager<IAtom[]> bondingAgentManager) {
        super();
        simulation = sim;
        ISpecies lastSpecies = sim.getSpecies(sim.getSpeciesCount() - 1);
        AtomType lastAtomType = lastSpecies.getAtomType(lastSpecies.getAtomTypeCount() - 1);
        colorMaps = new Color[lastAtomType.getIndex() + 1][4];
        this.bondingAgentManager = bondingAgentManager;
    }

    public Color getAtomColor(IAtom atom) {
        IAtom[] nbrs = bondingAgentManager.getAgent(atom);
        if (nbrs != null) {
            return colorMaps[atom.getType().getIndex()][nbrs.length];
        }
        // we weren't told how to deal with any atom type with this many bonds.
        return ColorScheme.DEFAULT_ATOM_COLOR;
    }

    /**
     * Sets atoms of the given type and number of bonds to be the given color.
     */
    public void setColor(AtomType type, int nBonds, Color color) {
        colorMaps[type.getIndex()][nBonds] = color;
    }
}
