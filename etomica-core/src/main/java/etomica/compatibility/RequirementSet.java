/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.compatibility;

import java.io.Serializable;
import java.util.Iterator;
import java.util.LinkedList;

public class RequirementSet extends Requirement implements Serializable
{
	RequirementSet() {}
	public RequirementSet add( Requirement req )
	{
		requirements.add( req );
		return this;
	}
	public boolean isSatisfied( FeatureSet featlist )
	{
		Iterator it = requirements.iterator();
		while ( it.hasNext() )
		{
			Requirement req = (Requirement) it.next();
			if ( !req.isSatisfied( featlist ) )
				return false;
		}
		return true;
	}
	private LinkedList requirements = new LinkedList();
};