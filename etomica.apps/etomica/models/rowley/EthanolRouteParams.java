/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.models.rowley;

/* 
 * Parameters to specify the approach route for dimer of ethanol molecules.
 * Please see Rowley et al (2006) for details.
 * 
 * K.R. Schadel 2008
 */

public class EthanolRouteParams {

	public static double[][] setEthanolParams(int route) {
		
		double[][] params;
				
		                              //   yaw      pitch      roll       theta      phi                     liftOff
		if (route == 1) {
			params = new double [][] { {   0.00}, { -73.49}, { 180.00}, {-180.00}, {-1.5156, 1.1837, 0.7967} };
		} else if (route == 2) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {-180.00}, {-0.071 ,-0.0296, 0.801 } };
		} else if (route == 3) {
			params = new double [][] { { 180.00}, { -71.62}, {  -0.04}, {   0.00}, {-0.8611, 1.8109,-0.5978} };
		} else if (route == 4) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {   0.00}, {-0.0025,-0.0291,-0.8015} };
		} else if (route == 5) {
			params = new double [][] { {   0.00}, {   0.00}, {   0.00}, {-180.00}, {-1.2419, 0.9066, 0.8086} };
		} else if (route == 6) {
			params = new double [][] { { 180.00}, {  36.76}, {   0.00}, {-180.00}, { 0.0263,-0.0306, 0.801 } };
		} else if (route == 7) {
			params = new double [][] { { 180.00}, {   0.00}, {-180.00}, {   0.00}, { 0.1915,-0.6075,-0.5862} };
		} else if (route == 8) {
			params = new double [][] { {   0.00}, {  71.63}, { 180.00}, {-180.00}, {-0.7462,-0.6323, 0.8008} };
		} else if (route == 9) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {-180.00}, {-1.2306,-0.4871, 0.7989} };
		} else if (route == 10) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {-180.00}, {-1.3661,-2.671 , 0.8373} };
		} else if (route == 11) {
			params = new double [][] { {   0.00}, {  18.48}, { 180.00}, {   0.00}, { 1.1532,-0.9933,-0.5905} };
		} else if (route == 12) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {-180.00}, { 1.4246,-1.2416,-0.3008} };
		} else if (route == 13) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {-180.00}, {-1.9507,-1.6644, 0.8139} };
		} else if (route == 14) {
			params = new double [][] { { 180.00}, { -18.45}, {   0.00}, {-180.00}, {-2.7077, 1.7952, 0.5955} };
		} else if (route == 15) {
			params = new double [][] { {   0.00}, {  42.97}, { 180.00}, {   0.00}, { 1.4663, 0.2158, 0.3255} } ;
		} else if (route == 16) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {-180.00}, { 0.7954, 0.4257, 0.8045} };
		} else if (route == 17) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {   0.00}, {-2.3068,-1.3416, 0.3212} };
		} else if (route == 18) {
			params = new double [][] { { 180.00}, {   0.00}, {-180.00}, {   0.00}, { 1.6928,-1.5437,-0.5876} };
		} else if (route == 19) {
			params = new double [][] { { 180.00}, {   0.00}, { 180.00}, {   0.00}, {18.377}                  };
		} else if (route == 20) {
			params = new double [][] { { 180.00}, {  36.75}, {   0.00}, {   0.00}, { 1.2828, 1.3308, 0.3132} };
		} else if (route == 21) {
			params = new double [][] { {   0.00}, {   1.06}, {   0.00}, {   0.00}, {-2.8989, 2.2807, 0.3306} };
		} else if (route == 22) {
			params = new double [][] { {   0.00}, { -71.62}, {   0.00}, {   0.00}, {-0.2772, 0.3586, 0.3179} };
		} else {
			throw new IllegalArgumentException("Must select route from [1,22].");
		}
		
		double pitch = params[1][0];
		
		if (route == 1) {
			
		} else if (route >=4 && route <= 5) {
			pitch = -pitch;
		} else if (route == 15) {
			// pitch = - (270 + pitch);
		}
		
		params[1][0] = pitch;
		
		return params;
	}
	
}
