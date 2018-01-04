/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.AtomType;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.AtomTypeAgentManager.AgentSource;
import etomica.atom.IAtom;
import etomica.simulation.Simulation;

import java.awt.*;
import java.util.HashMap;
import java.util.Map;

/**
 * Colors the atom according to the color given by its type field.
 *
 * @author David Kofke
 */
public class ColorSchemeByType extends ColorScheme {

    protected final Color[] moreDefaultColors = {ColorScheme.DEFAULT_ATOM_COLOR, Color.BLUE, Color.GREEN, Color.YELLOW, Color.ORANGE};
    private final Map<AtomType, Color> colorMap;
    protected int defaultColorsUsed = 0;

    public ColorSchemeByType() {
    	super();
        colorMap = new HashMap<>();
    }

    public void setColor(AtomType type, Color c) {
        colorMap.put(type, c);
    }

    public Color getAtomColor(IAtom a) {
        return getColor(a.getType());
    }

    public Color getColor(AtomType type) {
        return this.colorMap.computeIfAbsent(type, atomType -> {
            Color c;
            if (defaultColorsUsed < moreDefaultColors.length) {
                c = moreDefaultColors[defaultColorsUsed];
                defaultColorsUsed++;
            } else {
                c = defaultColor;
            }
            return c;
        });
    }
}
