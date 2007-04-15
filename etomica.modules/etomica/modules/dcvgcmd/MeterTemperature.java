package etomica.modules.dcvgcmd;

import etomica.EtomicaInfo;
import etomica.atom.AtomArrayList;
import etomica.atom.IAtomGroup;
import etomica.atom.SpeciesAgent;
import etomica.species.Species;
import etomica.units.Dimension;
import etomica.units.Temperature;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition
 */
public final class MeterTemperature extends etomica.data.meter.MeterTemperature {

	public MeterTemperature(Species species) {
		super();
		this.species = species;
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records temperature as given via kinetic energy");
		return info;
	}

	public double getDataAsScalar() {
		SpeciesAgent agent = phase.getAgent(species);
		AtomArrayList list = agent.getChildList();
		int size = list.size();
		int natoms = 0;
		if(size > 0) {
			natoms = size * ((IAtomGroup)list.get(0)).getChildList().size();
		}
		return (2. / ((phase.atomCount()- natoms) * phase.getBoundary().getDimensions().getD()))
				* meterKE.getDataAsScalar();
	}

	public Dimension getDimension() {
		return Temperature.DIMENSION;
	}

    private static final long serialVersionUID = 1L;
	private final Species species;
}