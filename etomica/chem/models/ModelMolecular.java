/*
 * Created on Jan 16, 2004
 */
package etomica.chem.models;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryHetero;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomPositionDefinitionSimple;
import etomica.atom.AtomTreeNodeFactory;
import etomica.atom.AtomTreeNodeGroupArray;
import etomica.atom.AtomTypeGroup;
import etomica.chem.Model;
import etomica.config.Conformation;
import etomica.simulation.Simulation;
import etomica.species.Species;
/**
 * Model for a molecule, which is formed from one or more sub-models, which
 * typically would describe the atoms in the molecule, or groups of atoms.
 */
public abstract class ModelMolecular extends Model {
	
	public final Conformation conformation;
	public final Model[] models;
	public final int[] count;
	
	public ModelMolecular(Conformation c, Model[] models) {
		this(c, models, defaultCount(models.length));//default is one of each model in the array
	}
	
	public ModelMolecular(Conformation c, Model[] models, int[] count) {
		super();
		if(models.length != count.length) throw new IllegalArgumentException("Inconsistent sizes of model and count array in ModelMolecular constructor");
		conformation = c;
		this.models = models;
		this.count = count;
		setDoNeighborIteration(true);		
	}
	
	//returns an array of 1's of length k
	private static int[] defaultCount(int k) {
		int[] c = new int[k];
		for(int i=0; i<k; i++) c[i] = 1;
		return c;
	}
	/**
	 * Invokes superclass method then sets neighbor-iteration flag for all
	 * submodels to false.
	 * @see etomica.chem.Model#setDoNeighborIteration(boolean)
	 */
	public void setDoNeighborIteration(boolean doNeighborIteration) {
		super.setDoNeighborIteration(doNeighborIteration);
		if(doNeighborIteration) {//set all sub-models to not iterate neighbors if this one is
			for(int i=0; i<models.length; i++) {
				models[i].setDoNeighborIteration(true);//need to set direct children to true to cascade a "false" setting to its children
				models[i].setDoNeighborIteration(false);//then set direct children to false too
			} 
		}
	}
	
	public AtomFactory makeAtomFactory(Simulation sim) {
		AtomTreeNodeFactory nodeFactory = AtomTreeNodeGroupArray.FACTORY;
//		switch(childCount) {
//			case 3: nodeFactory = AtomTreeNode3Site.FACTORY; break;
//			default: nodeFactory = AtomTreeNodeGroup.FACTORY; break;
//		}

		if(models.length == 1) {
			AtomFactoryHomo factory = new AtomFactoryHomo(sim.space, 
                    new AtomTypeGroup(new AtomPositionDefinitionSimple()), 
                    nodeFactory, count[0], conformation);
            factory.getType().setParentType(Species.makeAgentType(sim));
            factory.setChildFactory(models[0].makeAtomFactory(sim));
            return factory;
		}
		//makes array of child factories from models and counts
        int childCount = 0;
        for(int i=0; i<count.length; i++) childCount += count[i];//total number of child atoms in group
		AtomFactoryHetero factory = new AtomFactoryHetero(sim.space, 
                new AtomTypeGroup(new AtomPositionDefinitionSimple()), 
                nodeFactory, conformation);
        factory.getType().setParentType(Species.makeAgentType(sim));
        AtomFactory[] childFactories = new AtomFactory[childCount];
        for(int i=0,k=0; i<models.length; i++) {
            AtomFactory childFactory = models[i].makeAtomFactory(sim);
            for(int j=0; j<count[i]; j++) childFactories[k++] = childFactory;
        }
        return factory;
	}

}


