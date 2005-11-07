package etomica.modules.dcvgcmd;

import etomica.EtomicaInfo;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.atom.SpeciesAgent;
import etomica.species.Species;
import etomica.units.Dimension;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition
 */

/*
 * History of changes 7/03/02 (DAK) Changes to tie in with function of
 * kinetic-energy meter.
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
		AtomList list = ((AtomTreeNodeGroup)agent.node).childList;
		int size = list.size();
		int natoms = 0;
		if(size > 0) {
			natoms = size * ((AtomTreeNodeGroup)list.getFirst().node).childList.size();
		}
		return (2. / (double) ((phase.atomCount()- natoms) * phase.getBoundary().dimensions().D()))
				* meterKE.getDataAsScalar();
	}

	public Dimension getDimension() {
		return Dimension.TEMPERATURE;
	}

	private final Species species;
}