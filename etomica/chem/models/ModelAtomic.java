/*
 * Created on Jan 16, 2004
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package etomica.chem.models;
import etomica.AtomFactory;
import etomica.Simulation;
import etomica.Species;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomSequencerFactory;
import etomica.atom.AtomTypeSphere;
import etomica.chem.Electrostatic;
import etomica.chem.Element;
import etomica.chem.Model;

/**
 * @author zhaofang
 *
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
public abstract class ModelAtomic extends Model {

	private final Element element;
	private final Electrostatic electrostatic;
	
	public ModelAtomic() {
		this(new etomica.chem.elements.Undefined()); 
	}
	
	public ModelAtomic(Element element) {
		this(element, null);
	}
	
	public ModelAtomic(Element element, Electrostatic electrostatic) {
		super();
		this.element = element;
		this.electrostatic = electrostatic;
		setDoNeighborIteration(true);	
	}
	
	public AtomFactory makeAtomFactory(Simulation sim) {
        AtomSequencerFactory seqFactory = doNeighborIteration() ? sim.potentialMaster.sequencerFactory()
                 : AtomSequencerFactory.SIMPLE;
		return new AtomFactoryMono(sim.space,new AtomTypeSphere(Species.makeAgentType(sim)),seqFactory);
	}
	/**
	 * Returns the electrostatic.
	 * @return Electrostatic
	 */
	public Electrostatic getElectrostatic() {
		return electrostatic;
	}

	/**
	 * Returns the element.
	 * @return Element
	 */
	public Element getElement() {
		return element;
	}

}
