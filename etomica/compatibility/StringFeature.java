package etomica.compatibility;

import etomica.compatibility.Feature;

	public class StringFeature extends Feature
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