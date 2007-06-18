package etomica.paracetamol;

import java.lang.reflect.Constructor;

import etomica.EtomicaElement;
import etomica.simulation.ISimulation;
import etomica.simulation.Simulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

public class SpeciesParacetamol extends Species implements EtomicaElement {

	public SpeciesParacetamol(ISimulation sim) {
		super(new AtomFactoryParacetamol(sim));
	}
	
	public SpeciesSignature getSpeciesSignature(){
		Constructor constructor = null;
		try {
			constructor = this.getClass().getConstructor(new Class[] {Simulation.class});
		}
		catch(NoSuchMethodException e){
			System.err.println("Have NO CONSTRUCTOR");
		}
		return new SpeciesSignature(getName(), constructor,new Object[]{});
	}
	
	private static final long serialVersionUID = 1L;
}