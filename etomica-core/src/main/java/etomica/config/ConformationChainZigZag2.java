/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.config;

import etomica.api.IVector;
import etomica.space.Space;

public class ConformationChainZigZag2 extends ConformationChain {
	
    public ConformationChainZigZag2(Space space){
		super(space);
		
		v1 = space.makeVector();
		v2 = space.makeVector();
		isVectorOne = true;
	}
	public ConformationChainZigZag2(Space space, IVector vect1, IVector vect2){
		super(space);
		v1 = space.makeVector();
		v2 = space.makeVector();
		
		v1.E(vect1);
		v2.E(vect2);
		
		isVectorOne = true;
	}
	
	/* (non-Javadoc)
	 * @see etomica.ConformationChain#reset()
	 */
	protected void reset() {
		isVectorOne = true;
	}

	/* (non-Javadoc)
	 * @see etomica.ConformationChain#nextVector()
	 */
	protected IVector nextVector() {
		if(isVectorOne){
			isVectorOne = false;
			return v1;
		}
		isVectorOne = true;
		return v2;
	}
	
	public IVector getFirstVector() {
	    return v1;
	}
	
	public IVector getSecondVector() {
	    return v2;
	}

    private static final long serialVersionUID = 1L;
	protected IVector v1;
	protected IVector v2;
	boolean isVectorOne;
   }


