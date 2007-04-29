package etomica.compatibility;

import java.io.Serializable;

public abstract class Requirement implements Serializable
{
	public abstract boolean isSatisfied( FeatureSet featlist );
}