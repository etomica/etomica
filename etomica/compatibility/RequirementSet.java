package etomica.compatibility;

import java.io.Serializable;
import java.util.Iterator;
import java.util.LinkedList;

import etomica.compatibility.Requirement;

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