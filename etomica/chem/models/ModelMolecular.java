/*
 * Created on Jan 16, 2004
 */
package etomica.chem.models;
import etomica.AtomFactory;
import etomica.AtomTreeNode;
import etomica.AtomTreeNodeFactory;
import etomica.Configuration;
import etomica.Space;
import etomica.atom.AtomFactoryHetero;
import etomica.atom.AtomFactoryHomo;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTreeNodeGroupArray;
import etomica.chem.Model;
/**
 * Model for a molecule, which is formed from one or more sub-models, which
 * typically would describe the atoms in the molecule, or groups of atoms.
 */
public abstract class ModelMolecular extends Model {
	
	public final Configuration configuration;
	public final Model[] models;
	public final int[] count;
	
	public ModelMolecular(Configuration c, Model[] models) {
		this(c, models, defaultCount(models.length));//default is one of each model in the array
	}
	
	public ModelMolecular(Configuration c, Model[] models, int[] count) {
		super();
		if(models.length != count.length) throw new IllegalArgumentException("Inconsistent sizes of model and count array in ModelMolecular constructor");
		configuration = c;
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
	
	public AtomFactory makeAtomFactory(Space space) {
		AtomSequencerFactory seqFactory = doNeighborIteration() ? sim.iteratorFactory.neighborSequencerFactory()
																 : sim.iteratorFactory.simpleSequencerFactory();
		int childCount = 0;
		for(int i=0; i<count.length; i++) childCount += count[i];//total number of child atoms in group
		AtomTreeNodeFactory nodeFactory = AtomTreeNodeGroupArray.FACTORY;
//		switch(childCount) {
//			case 3: nodeFactory = AtomTreeNode3Site.FACTORY; break;
//			default: nodeFactory = AtomTreeNodeGroup.FACTORY; break;
//		}

		if(models.length == 1) {
			AtomFactory childFactory = models[0].makeAtomFactory(space);
			return new AtomFactoryHomo(space, seqFactory, nodeFactory, 
										childFactory, count[0], configuration);
		}
		//makes array of child factories from models and counts
		AtomFactory[] childFactories = new AtomFactory[childCount];
		int k = 0;
		for(int i=0; i<models.length; i++) {
			AtomFactory childFactory = models[i].makeAtomFactory(space);
			for(int j=0; j<count[i]; j++) childFactories[k++] = childFactory;
		}
		return new AtomFactoryHetero(space, seqFactory, nodeFactory, childFactories, configuration);
	}

}


