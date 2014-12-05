/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

/*
 * Created on Jan 30, 2004
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package etomica.chem.electrostatics;

/**
 * @author zhaofang
 *
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
public class Charge implements Electrostatic, java.io.Serializable {
	double q;
	
	/**
	 * 
	 */
	public Charge() {
		this(0.0);
	}

	public Charge(double q) {
		this.q = q;	
	}
	/**
	 * @return
	 */
	public double getQ() {
		return q;
	}

	/**
	 * @param d
	 */
	public void setQ(double d) {
		q = d;
	}

}
