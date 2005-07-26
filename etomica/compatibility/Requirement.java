package etomica.compatibility;

import java.io.Serializable;

import etomica.compatibility.FeatureSet;

public abstract class Requirement implements Serializable
{
	public abstract boolean isSatisfied( FeatureSet featlist );
}