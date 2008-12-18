package etomica.graphics;
import java.awt.Color;

import etomica.api.IAtomLeaf;
import etomica.api.IBox;
import etomica.api.IMolecule;
import etomica.api.IRandom;
import etomica.api.ISimulation;
import etomica.atom.MoleculeAgentManager;
import etomica.atom.MoleculeAgentManager.MoleculeAgentSource;

public class ColorSchemeRandomByMolecule extends ColorScheme implements MoleculeAgentSource {
    
    public ColorSchemeRandomByMolecule(ISimulation sim, IBox box, IRandom random) {
    	super(sim);
        this.random = random;
        agentManager = new MoleculeAgentManager(sim, box, this);
    }
    
    public Color getAtomColor(IAtomLeaf a) {
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