/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.compatibility;

import java.io.Serializable;

public class ClassRequirement extends Requirement implements Serializable
{
	ClassRequirement( Class aclass )
	{
		comparison = new StringFeature( "CLASS", aclass.getName() );
	}
	public boolean isSatisfied( FeatureSet featlist )
	{
		if ( comparison==null )
			return false;
		int result = featlist.get( "CLASS" ).compareTo( comparison );
		return result==0; 
	}
	Feature comparison = null;
}
