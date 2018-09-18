/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;
import java.awt.Color;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomLeafAgentManager.AgentSource;

/**
 * Parent class for color schemes that are best implemented by attaching colors
 * to all atoms at once, rather than determining the color of each as it is drawn.
 * The colorAllAtoms method is called by the display if it determines that the
 * ColorScheme is a subclass of this one.
 */
public abstract class ColorSchemeCollectiveAgent extends ColorScheme implements AgentSource<Color>, ColorSchemeCollective {
    
    protected AtomLeafAgentManager<Color> agentManager;
    
    public ColorSchemeCollectiveAgent(Box box) {
    	super();
        agentManager = new AtomLeafAgentManager<Color>(this, box);
    }
    
    public abstract void colorAllAtoms();
    
    public Color getAtomColor(IAtom a) {
        return agentManager.getAgent(a);
    }
   
    public Class getAgentClass() {
        return Color.class;
    }
    
    public Color makeAgent(IAtom a, Box agentBox) {
        return null;
    }
    
    public void releaseAgent(Color agent, IAtom atom, Box agentBox) {}
    
}
