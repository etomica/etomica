/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

import etomica.atom.IAtom;
import etomica.box.Box;
import etomica.molecule.IMolecule;
import etomica.molecule.MoleculeAgentManager;
import etomica.molecule.MoleculeAgentManager.MoleculeAgentSource;
import etomica.species.SpeciesManager;
import etomica.util.random.IRandom;

import java.awt.*;

public class ColorSchemeRandomByMolecule extends ColorScheme implements MoleculeAgentSource<Color> {

    public ColorSchemeRandomByMolecule(SpeciesManager sm, Box box, IRandom random) {
        super();
        this.random = random;
        agentManager = new MoleculeAgentManager<>(sm, box, this);
    }

    public Color getAtomColor(IAtom a) {
        return agentManager.getAgent(a.getParentGroup());
    }

    public Color makeAgent(IMolecule a) {
        return new Color((float) random.nextDouble(), (float) random.nextDouble(), (float) random.nextDouble());
    }

    public void releaseAgent(Color agent, IMolecule atom) {
    }

    private final MoleculeAgentManager<Color> agentManager;
    private final IRandom random;
}
