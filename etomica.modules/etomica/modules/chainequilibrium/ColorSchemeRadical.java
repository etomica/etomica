package etomica.modules.chainequilibrium;

import java.awt.Color;

import etomica.api.IAtom;
import etomica.api.IAtomTypeLeaf;
import etomica.api.ISimulation;
import etomica.atom.AtomLeafAgentManager;
import etomica.atom.AtomTypeAgentManager;
import etomica.graphics.ColorSchemeByType;

public class ColorSchemeRadical extends ColorSchemeByType {

    public ColorSchemeRadical(ISimulation sim, AtomLeafAgentManager agentManager) {
        super(sim);
        this.agentManager = agentManager;
        radicalColorMap = new AtomTypeAgentManager(this, sim.getSpeciesManager(),
                sim.getEventManager(), false);
        fullColorMap = new AtomTypeAgentManager(this, sim.getSpeciesManager(),
                sim.getEventManager(), false);
    }
    
    public Color getAtomColor(IAtom atom) {
        IAtom[] nbrs = (IAtom[])agentManager.getAgent(atom);
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
            return super.getAtomColor(atom);
        }
        return (Color)radicalColorMap.getAgent(atom.getType());
    }

    public void setFreeRadicalColor(IAtomTypeLeaf type, Color color) {
        radicalColorMap.setAgent(type, color);
    }
    
    public void setFullColor(IAtomTypeLeaf type, Color color) {
        fullColorMap.setAgent(type, color);
    }

    protected final AtomLeafAgentManager agentManager;
    protected final AtomTypeAgentManager radicalColorMap, fullColorMap;
}
