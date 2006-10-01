/*
 * Created on Jan 16, 2004
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package etomica.chem.models;
import etomica.atom.AtomFactory;
import etomica.atom.AtomFactoryMono;
import etomica.atom.AtomTypeSphere;
import etomica.chem.electrostatics.Electrostatic;
import etomica.chem.elements.Element;
import etomica.simulation.Simulation;
import etomica.space.CoordinateFactorySphere;
import etomica.species.Species;

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
        this(null);
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
        AtomFactoryMono factory = new AtomFactoryMono(new CoordinateFactorySphere(sim),new AtomTypeSphere(sim));
        factory.getType().setParentType(Species.makeAgentType(sim));
        return factory;
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
