package etomica.compatibility;

import java.io.Serializable;


/** Tentative property-based infrastructure for a general property-based compatibility algorithm */
public abstract class Feature implements Serializable
{
	public static final int LESS_THAN = -1;
	public static final int IS_EQUAL = 0;
	public static final int GREATER_THAN = 1;
	public static final int IS_EMPTY = 2;
	public static final int INCOMPATIBLE_TYPES = 3;
	public static final int CONTAINS = 4;
	public Feature( String aname ) { name=aname; }
	public abstract int compareTo( Feature feat );
	public abstract boolean compareTo( int operator, Feature feat );
	public String name;
};