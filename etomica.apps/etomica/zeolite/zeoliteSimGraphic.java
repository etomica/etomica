package etomica.zeolite;

import etomica.graphics.SimulationGraphic;
import etomica.space.Space;

public class zeoliteSimGraphic extends SimulationGraphic {

	public zeoliteSimGraphic(ZeoliteSimulation sim, Space space) {
		this(sim, space, "");
	}

	public zeoliteSimGraphic(ZeoliteSimulation sim, Space space, String appName){
		super(sim, appName, space);
		DeviceButtonSingle tSwitch = new DeviceButtonSingle(sim.getController());
		ZeoliteSimStart action = new ZeoliteSimStart(sim,this);
		tSwitch.setAction(action);
		tSwitch.setLabel("Start Simulation");
		super.add(tSwitch);
		
	}
	
}
