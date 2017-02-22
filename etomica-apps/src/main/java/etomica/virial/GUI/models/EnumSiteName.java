/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.GUI.models;

public enum EnumSiteName {
	
	C("Carbon"),
	O("Oxygen"),
	H("Hydrogen"),
	CH3("Methyl"),
	CH2("Methylene"),
	CH("Methyne"),
	OH("Hydroxyl"),
	aC("AlphaCarbon"),
	aH("AlphaHydrogen"),
	X("Satellite"),
	LJ("Leaf");
	
	
	
	
	private String Site;
	
	public String getSite() {
		return Site;
	}

	public void setSite(String site) {
		Site = site;
	}

	EnumSiteName(String sitename){
		Site = sitename;
	}
	

}
