/*
 * Created on Jan 30, 2004
 *
 * To change the template for this generated file go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
package etomica.chem.electrostatics;
import etomica.chem.*;

/**
 * @author zhaofang
 *
 * To change the template for this generated type comment go to
 * Window>Preferences>Java>Code Generation>Code and Comments
 */
public class Charge implements Electrostatic {
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
