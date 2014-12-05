/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.virial.cluster;

import etomica.virial.ClusterBonds;
import etomica.virial.MayerFunction;

/**
 * @author kofke
 *
 * A cluster formed with a bond between every pair.
 */
public class Full extends ClusterBonds {

	/**
	 * Constructor for a fully joined cluster
	 * @param n number of points in ring
	 */
	public Full(int n) {
		super(n, new int[][][] {Standard.full(n)});
	}
}
