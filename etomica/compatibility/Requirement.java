package etomica.compatibility;

import etomica.compatibility.FeatureSet;

public abstract class Requirement
{
	public abstract boolean isSatisfied( FeatureSet featlist );
}