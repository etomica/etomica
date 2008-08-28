package etomica.zeolite;

import etomica.graphics.SimulationGraphic;
import etomica.space.ISpace;

public class zeoliteSimGraphic extends SimulationGraphic {

	public zeoliteSimGraphic(ZeoliteSimulation sim, ISpace space) {
		this(sim, space, "");
	}

	public zeoliteSimGraphic(ZeoliteSimulation sim, ISpace space, String appName){
		super(sim, appName, space, sim.getController());
		DeviceButtonSingle tSwitch = new DeviceButtonSingle(sim.getController());
		ZeoliteSimStart action = new ZeoliteSimStart(sim,this);
		tSwitch.setAction(action);
		tSwitch.setLabel("Start Simulation");
		super.add(tSwitch);
		
	}
	
}
