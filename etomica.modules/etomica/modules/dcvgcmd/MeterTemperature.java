package etomica.modules.dcvgcmd;

import etomica.EtomicaElement;
import etomica.EtomicaInfo;
import etomica.Phase;
import etomica.Species;
import etomica.SpeciesAgent;
import etomica.atom.AtomList;
import etomica.atom.AtomTreeNodeGroup;
import etomica.data.meter.MeterKineticEnergy;
import etomica.data.meter.MeterScalar;
import etomica.units.Dimension;

/**
 * Meter for measurement of the temperature based on kinetic-energy
 * equipartition
 */

/*
 * History of changes 7/03/02 (DAK) Changes to tie in with function of
 * kinetic-energy meter.
 */

public final class MeterTemperature extends MeterScalar implements
		EtomicaElement {

	public MeterTemperature(Species species) {
		super();
		this.species = species;
		setLabel("Temperature");
		meterKE = new MeterKineticEnergy();
	}

	public static EtomicaInfo getEtomicaInfo() {
		EtomicaInfo info = new EtomicaInfo(
				"Records temperature as given via kinetic energy");
		return info;
	}

	public double getDataAsScalar(Phase phase) {
		SpeciesAgent agent = phase.getAgent(species);
		AtomList list = ((AtomTreeNodeGroup)agent.node).childList;
		int size = list.size();
		int natoms = 0;
		if(size > 0){
			natoms = size * ((AtomTreeNodeGroup)list.getFirst().node).childList.size();
		}
		return (2. / (double) ((phase.atomCount()- natoms) * phase.boundary().dimensions().D()))
				* meterKE.getDataAsScalar(phase);
	}

	public Dimension getDimension() {
		return Dimension.TEMPERATURE;
	}

	private final MeterKineticEnergy meterKE;
	private final Species species;
}