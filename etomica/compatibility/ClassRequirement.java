package etomica.compatibility;

import java.io.Serializable;

import etomica.compatibility.Requirement;

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
