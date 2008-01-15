package etomica.paracetamol;

import java.lang.reflect.Constructor;

import etomica.simulation.ISimulation;
import etomica.species.Species;
import etomica.species.SpeciesSignature;

public class SpeciesParacetamol extends Species {

	public SpeciesParacetamol(ISimulation sim) {
		super();
		setMoleculeFactory(new AtomFactoryParacetamol(sim, this));
	}
	
	public SpeciesSignature getSpeciesSignature(){
		Constructor constructor = null;
		try {
			constructor = this.getClass().getConstructor(new Class[] {ISimulation.class});
		}
		catch(NoSuchMethodException e){
			System.err.println("Have NO CONSTRUCTOR");
		}
		return new SpeciesSignature(constructor,new Object[]{});
	}
	
	private static final long serialVersionUID = 1L;
}