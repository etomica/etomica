/*
 * Created on Jan 16, 2004
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package etomica.chem.models;
import etomica.Space;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
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
	
	public AtomFactory makeAtomFactory(Space space) {
		AtomSequencer.Factory seqFactory = doNeighborIteration() ? sim.iteratorFactory.neighborSequencerFactory()
																 : sim.iteratorFactory.simpleSequencerFactory();
		return new AtomFactoryMono(space,seqFactory);
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
