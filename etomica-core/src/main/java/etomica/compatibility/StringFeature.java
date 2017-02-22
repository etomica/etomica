/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.compatibility;

import java.io.Serializable;

	public class StringFeature extends Feature implements Serializable
	{
		public StringFeature( String aname, String avalue ) { super(aname); value=avalue; }
		public int compareTo( Feature feat )
		{
			if ( feat==null)
				return IS_EMPTY;
			if ( feat instanceof StringFeature )
				return ((StringFeature)feat).value.compareTo( value );
			return INCOMPATIBLE_TYPES;		
		}
		public boolean compareTo( int operator, Feature feat )
		{
			if ( operator == Feature.CONTAINS )
				return ((StringFeature) feat).value.indexOf( value )!=-1;
			return compareTo( feat )==operator;
		}

		public String value;
	}