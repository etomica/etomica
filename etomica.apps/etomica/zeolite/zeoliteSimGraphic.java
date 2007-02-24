package etomica.zeolite;

import etomica.graphics.DisplayPhase;
import etomica.graphics.SimulationGraphic;

public class zeoliteSimGraphic extends SimulationGraphic{
	public DisplayPhase display;
	public zeoliteSimGraphic(ZeoliteSimulation sim){
		super(sim);
		DeviceButtonSingle tSwitch = new DeviceButtonSingle(sim.getController());
		ZeoliteSimStart action = new ZeoliteSimStart(sim,this);
		tSwitch.setAction(action);
		tSwitch.setLabel("Start Simulation");
		super.add(tSwitch);
		
		display = super.getDisplayPhase(sim.phase);
		/*
		DeviceButton picture = new DeviceButton(sim.getController());
		TakePictures pic = new TakePictures(this);
		picture.setAction(pic);
		picture.setLabel("Picture Time");
		super.add(picture);
		*/
	}
	
}
