package etomica.zeolite;

import etomica.graphics.SimulationGraphic;

public class zeoliteSimGraphic extends SimulationGraphic {

	public zeoliteSimGraphic(ZeoliteSimulation sim) {
		this(sim, "");
	}

	public zeoliteSimGraphic(ZeoliteSimulation sim, String appName){
		super(sim, appName);
		DeviceButtonSingle tSwitch = new DeviceButtonSingle(sim.getController());
		ZeoliteSimStart action = new ZeoliteSimStart(sim,this);
		tSwitch.setAction(action);
		tSwitch.setLabel("Start Simulation");
		super.add(tSwitch);
		
	}
	
}
