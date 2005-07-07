package etomica.compatibility;

import etomica.EtomicaInfo;
import etomica.Potential;
import etomica.compatibility.Feature;
import etomica.compatibility.Requirement;
import etomica.potential.P1BondedHardSpheres;
import etomica.simulations.HSMD3D;
import etomica.space2d.Space2D;

import java.util.HashMap;


/** Container for a class features. EtomicaInfo has a method getFeatures() that returns a FeatureSet.
 * Generic algorithms may use the Requirements object interface to test an object's feature for compatibility.
 * @author Henrique
 */
public final class FeatureSet
{
	public FeatureSet() {}
	public FeatureSet add( Feature feat ) { list.put( feat.name, feat ); return this; }
	public FeatureSet add( String name, String value ) { list.put( name,new StringFeature(name,value) ); return this; }
	public FeatureSet add( String name, double value ) { list.put( name,new NumericFeature(name,value) ); return this; }
	public Feature get( String name ) { return (Feature) list.get( name ); }
	public final boolean satisfies( Requirement reqs )	{ return reqs.isSatisfied( this );	}
	protected HashMap list = new HashMap();
	
	static public void main( String[] args )
	{
		//Potential potential = new P1BondedHardSpheres( new Space2D() );
		RequirementSet req = new RequirementSet()
		.add( new FeatureRequirement( "SPACEDIM", Feature.GREATER_THAN, 1 ) )
		.add( new FeatureRequirement( "NBODIES", Feature.IS_EQUAL, 1 ) );
		
		FeatureSet fss = EtomicaInfo.getInfo( P1BondedHardSpheres.class ).getFeatures();
		if ( fss.satisfies( req ) )
			System.out.println( "Potential  SATISFIED requiments." );
		else
			System.out.println( "Potential  DID NOT SATISFY requiments." );
	}
}
