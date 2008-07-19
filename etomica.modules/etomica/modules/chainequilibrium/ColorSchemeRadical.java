package etomica.modules.chainequilibrium;

import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IAtomLeaf;
import etomica.api.IAtomSet;
import etomica.api.IAtomTypeLeaf;
import etomica.api.IBox;
import etomica.api.ISimulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeAgentManager;
import etomica.graphics.ColorScheme;
import etomica.graphics.ColorSchemeByType;
import etomica.graphics.ColorSchemeCollective;

public class ColorSchemeRadical extends ColorSchemeByType implements ColorSchemeCollective, AtomLeafAgentManager.AgentSource {

    public ColorSchemeRadical(ISimulation sim, AtomLeafAgentManager agentManager) {
        super(sim);
        this.agentManager = agentManager;
        radicalColorMap = new AtomTypeAgentManager(this, sim.getSpeciesManager(),
                sim.getEventManager(), false);
        fullColorMap = new AtomTypeAgentManager(this, sim.getSpeciesManager(),
                sim.getEventManager(), false);
        greys = new Color[20];
        for (int i=0; i<20; i++) {
            greys[i] = new Color(90+6*i,90+6*i,90+6*i);
        }
    }

    public Color getAtomColor(IAtom atom) {
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(atom);
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
            int chainNumber = ((LengthAgent)chainLengthManager.getAgent(atom)).chainNumber % greys.length;
            return greys[chainNumber];
        }
        return (Color)radicalColorMap.getAgent(atom.getType());
    }

    public void setFreeRadicalColor(IAtomTypeLeaf type, Color color) {
        radicalColorMap.setAgent(type, color);
    }
    
    public void setFullColor(IAtomTypeLeaf type, Color color) {
        fullColorMap.setAgent(type, color);
    }
    
    public Class getAgentClass() {
        return LengthAgent.class;
    }

    public Object makeAgent(IAtomLeaf a) {
        return new LengthAgent();
    }

    public void releaseAgent(Object agent, IAtomLeaf atom) {}

    public void colorAllAtoms() {
        // untag all the Atoms
        IAtomSet leafList = box.getLeafList();
        int nLeaf = leafList.getAtomCount();
        for (int i=0; i<nLeaf; i++) {
            ((LengthAgent)chainLengthManager.getAgent(leafList.getAtom(i))).chainNumber = -1;
        }

        int chainNumber = 0;
        for (int i=0; i<nLeaf; i++) {
            IAtomLeaf a = (IAtomLeaf)leafList.getAtom(i);
            // if an Atom has a chain length, it was already counted as part of 
            // another chain
            if (((LengthAgent)chainLengthManager.getAgent(a)).chainNumber > 0) continue;

            recursiveTag(a, chainNumber);
            chainNumber++;
        }
    }

    protected void recursiveTag(IAtomLeaf a, int chainNumber) {
        ((LengthAgent)chainLengthManager.getAgent(a)).chainNumber = chainNumber;

        IAtomLeaf[] nbrs = (IAtomLeaf[])agentManager.getAgent(a);

        // count all the bonded partners
        for(int i=0; i<nbrs.length; i++) {
            if(nbrs[i] == null) continue;
            if(((LengthAgent)chainLengthManager.getAgent(nbrs[i])).chainNumber > -1) {
                // this Atom was already counted as being within this chain
                // so skip it
                continue;
            }
            // count this Atom and all of its bonded partners
            recursiveTag(nbrs[i], chainNumber);
        }
        return;
    }

    public IBox getBox() {
        return box;
    }

    public void setBox(IBox box) {
        this.box = box;
        if (chainLengthManager != null) {
            // allow old agentManager to de-register itself as a BoxListener
            chainLengthManager.dispose();
        }
        chainLengthManager = new AtomLeafAgentManager(this,box);
    }

    public static class LengthAgent {
        public int chainNumber;
    }

    protected IBox box;
    protected AtomLeafAgentManager chainLengthManager;
    protected final AtomLeafAgentManager agentManager;
    protected final AtomTypeAgentManager radicalColorMap, fullColorMap;
    protected final Color[] greys;
}
