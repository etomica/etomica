/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.zeolite;
import etomica.space3d.Space3D;


public class Zeolite {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		// Converts the simulation data into something I can use
		String inputFile = "32_1000000_0.00611_2000_WCA";
		//Cut data file down to multiple parts
		int cut = 5;
		DataCutter data = new DataCutter(inputFile,cut);
		
		for(int i=0;i<cut;i++){
			String file = inputFile+"_"+i;
			Converter(file);
			System.out.println("File "+i+" converted");
		}
		
		//Converter(inputFile,32,2);     
        System.out.println("Finished");
	}
	public static void Converter(String inputFile) {
		// TODO Auto-generated method stub
		String outputFile = inputFile+"_Result.txt";
		MSDProcessor proc = new MSDProcessor(Space3D.getInstance(),inputFile,outputFile);
		
		//proc.setDeltaTmax(1);
		proc.fillArrays();
		System.out.println("Converter done");
	}

}
