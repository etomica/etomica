/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

import etomica.virial.simulations.VirialRowleyAlcohol;

public class MassProductionVirialRowleyAlcohol {
	
	public static void main(String[] args) {

		int seed;
		
		if (args.length == 1) {
			seed = Integer.parseInt(args[0]);
		} else {
			seed = 0;
		}
		
		int temp = 700; // Kelvin
		
		for (int i=0; i<20; i++) {
			
			seed = seed + 1;
		
			//VirialRowleyAlcohol sim = new VirialRowleyAlcohol();
			String[] string = new String[4];
			
			//string = new String[] {"2", ""+temp, "10000"};
			string[0] = "2";
			string[1] = Integer.toString(temp);
			string[2] = "10000";
			string[3] = Long.toString(seed);
			
			VirialRowleyAlcohol.main(string);
		
		}
		
	}

}
