/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.compatibility;

import java.io.Serializable;


public class FeatureRequirement extends Requirement implements Serializable
{
	FeatureRequirement( int desired_result, Feature value )
	{
		expected_result = desired_result;
		reference = value;
	}
	FeatureRequirement( String name, int desired_result, double value )
	{
		expected_result = desired_result;
		reference = new NumericFeature( name, value );
	}
	FeatureRequirement( String name, int desired_result, String value )
	{
		expected_result = desired_result;
		reference = new StringFeature( name, value );
	}
	public boolean isSatisfied( FeatureSet featlist )
	{
		return reference.compareTo( expected_result, featlist.get( reference.name ) );
	}
	protected int expected_result;
	protected Feature reference;
}
