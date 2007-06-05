package etomica.compatibility;

import java.io.Serializable;
import java.util.HashMap;


/** Container for a class features. EtomicaInfo has a method getFeatures() that returns a FeatureSet.
 * Generic algorithms may use the Requirements object interface to test an object's feature for compatibility.
 * @author Henrique
 */
public final class FeatureSet implements Serializable
{
	public FeatureSet() {}
	public FeatureSet add( Feature feat ) { list.put( feat.name, feat ); return this; }
	public FeatureSet add( String name, String value ) { list.put( name,new StringFeature(name,value) ); return this; }
	public FeatureSet add( String name, double value ) { list.put( name,new NumericFeature(name,value) ); return this; }
	public Feature get( String name ) { return (Feature) list.get( name ); }
	public final boolean satisfies( Requirement reqs )	{ return reqs.isSatisfied( this );	}
	protected HashMap list = new HashMap();
	
}
