/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.modules.chainequilibrium;

import etomica.atom.*;
import etomica.box.Box;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.ColorSchemeCollective;
import etomica.simulation.Simulation;

import java.awt.*;

public class ColorSchemeRadical extends ColorSchemeByType implements ColorSchemeCollective, AtomLeafAgentManager.AgentSource<ColorSchemeRadical.LengthAgent> {

    protected final AtomLeafAgentManager<IAtom[]> agentManager;
    protected final AtomTypeAgentManager radicalColorMap, fullColorMap;
    protected final Color[] greys;
    protected Box box;
    protected AtomLeafAgentManager<LengthAgent> chainLengthManager;

    public ColorSchemeRadical(Simulation sim, AtomLeafAgentManager<IAtom[]> agentManager) {
        super(sim);
        this.agentManager = agentManager;
        radicalColorMap = new AtomTypeAgentManager(this, sim);
        fullColorMap = new AtomTypeAgentManager(this, sim);
        greys = new Color[20];
        for (int i=0; i<20; i++) {
            greys[i] = new Color(90+6*i,90+6*i,90+6*i);
        }
    }

    public Color getAtomColor(IAtom atom) {
        IAtom[] nbrs = agentManager.getAgent(atom);
        if (nbrs == null) {
            return ColorScheme.DEFAULT_ATOM_COLOR;
        }
        for(int i=0; i < nbrs.length-1; ++i){
            if (nbrs[i] == null) {
                return super.getAtomColor(atom);
            }
        }
        if (nbrs[nbrs.length-1] != null) {
            Color color = (Color)fullColorMap.getAgent(atom.getType());
            if (color != null) {
                return color;
            }
            int chainNumber = chainLengthManager.getAgent(atom).chainNumber % greys.length;
            return greys[chainNumber];
        }
        return (Color)radicalColorMap.getAgent(atom.getType());
    }

    public void setFreeRadicalColor(AtomType type, Color color) {
        radicalColorMap.setAgent(type, color);
    }

    public void setFullColor(AtomType type, Color color) {
        fullColorMap.setAgent(type, color);
    }

    public LengthAgent makeAgent(IAtom a, Box agentBox) {
        return new LengthAgent();
    }

    public void releaseAgent(LengthAgent agent, IAtom atom, Box agentBox) {}

    public void colorAllAtoms() {
        // untag all the Atoms
        IAtomList leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            chainLengthManager.getAgent(leafList.getAtom(i)).chainNumber = -1;
        }

        int chainNumber = 0;
        for (int i=0; i<nLeaf; i++) {
            IAtom a = leafList.getAtom(i);
            // if an Atom has a chain length, it was already counted as part of
            // another chain
            if (chainLengthManager.getAgent(a).chainNumber > 0) continue;

            recursiveTag(a, chainNumber);
            chainNumber++;
        }
    }

    protected void recursiveTag(IAtom a, int chainNumber) {
        chainLengthManager.getAgent(a).chainNumber = chainNumber;

        IAtom[] nbrs = agentManager.getAgent(a);

        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(chainLengthManager.getAgent(nbrs[i]).chainNumber > -1) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            recursiveTag(nbrs[i], chainNumber);
        }
        return;
    }

    public Box getBox() {
        return box;
    }

    public void setBox(Box box) {
        this.box = box;
        if (chainLengthManager != null) {
            // allow old agentManager to de-register itself as a BoxListener
            chainLengthManager.dispose();
        }
        chainLengthManager = new AtomLeafAgentManager<LengthAgent>(this,box,LengthAgent.class);
    }

    public static class LengthAgent {
        public int chainNumber;
    }
}
