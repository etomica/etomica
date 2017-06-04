/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import java.awt.Color;

import etomica.atom.IAtom;
import etomica.atom.IAtomType;
import etomica.simulation.Simulation;
import etomica.atom.AtomTypeAgentManager;
import etomica.atom.AtomTypeAgentManager.AgentSource;

/**
 * Colors the atom according to the color given by its type field.
 *
 * @author David Kofke
 */
public class ColorSchemeByType extends ColorScheme implements AgentSource {
    
    public ColorSchemeByType(Simulation sim) {
    	super();
        colorMap = new AtomTypeAgentManager(this, sim);
    }

    public Object makeAgent(IAtomType atom) {
    	return null;
    }

    public void releaseAgent(Object obj, IAtomType atom) {
    }

    public Class getSpeciesAgentClass() {
    	return Color.class;
    }
    
    public void setColor(IAtomType type, Color c) {
    	colorMap.setAgent(type, c);
    }
    
    public Color getAtomColor(IAtom a) {
        return getColor(a.getType());
    }
    
    public Color getColor(IAtomType type) {
        Color color = (Color)colorMap.getAgent(type);
        if (color == null) {
            if (defaultColorsUsed < moreDefaultColors.length) {
                color = moreDefaultColors[defaultColorsUsed];
                defaultColorsUsed++;
                setColor(type, color);
            }
            else {
                color = defaultColor;
                setColor(type, color);
            }
        }
        return color;
    }
    
    private final AtomTypeAgentManager colorMap;
    protected final Color[] moreDefaultColors = new Color[]{ColorScheme.DEFAULT_ATOM_COLOR, Color.BLUE, Color.GREEN, Color.YELLOW, Color.ORANGE};
    protected int defaultColorsUsed = 0;
}
