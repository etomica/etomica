/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;

import etomica.graphics.SimulationGraphic;
import etomica.space.Space;

public class zeoliteSimGraphic extends SimulationGraphic {

	public zeoliteSimGraphic(ZeoliteSimulation sim, Space space) {
		this(sim, space, "");
	}

	public zeoliteSimGraphic(ZeoliteSimulation sim, Space space, String appName){
		super(sim, appName);
		DeviceButtonSingle tSwitch = new DeviceButtonSingle(sim.getController());
		ZeoliteSimStart action = new ZeoliteSimStart(sim,this);
		tSwitch.setAction(action);
		tSwitch.setLabel("Start Simulation");
		super.add(tSwitch);
		
	}
	
}
