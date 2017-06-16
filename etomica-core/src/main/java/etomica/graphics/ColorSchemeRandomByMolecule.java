/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.simulation.Simulation;
import etomica.util.random.IRandom;

import java.awt.*;

public class ColorSchemeRandomByMolecule extends ColorScheme implements MoleculeAgentSource {
    
    public ColorSchemeRandomByMolecule(Simulation sim, Box box, IRandom random) {
    	super();
        this.random = random;
        agentManager = new MoleculeAgentManager(sim, box, this);
    }
    
    public Color getAtomColor(IAtom a) {
        return (Color)agentManager.getAgent(a.getParentGroup());
    }
   
    public Class getMoleculeAgentClass() {
        return Color.class;
    }
    
    public Object makeAgent(IMolecule a) {
        return new Color((float)random.nextDouble(),(float)random.nextDouble(),(float)random.nextDouble());
    }
    
    public void releaseAgent(Object agent, IMolecule atom) {}

    private static final long serialVersionUID = 2L;
    private final MoleculeAgentManager agentManager;
    private final IRandom random;
}
