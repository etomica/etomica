/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

public class DeviceBoxValueChangedEvent extends java.awt.AWTEvent {
	
	private int newValue;

	public DeviceBoxValueChangedEvent(DeviceBox source, int value) {
		super(source, RESERVED_ID_MAX+1);
		newValue = value;	
	} // end public DeviceBoxValueChangedEvent

	/**
	 * 
	 * @return the value of DeviceBox' text box.
	 */
	public int getValue() {
		return newValue;
	} // end getValue

	/**
	 *  Set the value of a DeviceBoxValueChangedEvent (should be the value of the text box)
	 * @param value
	 */
	public void setValue(int value) {
		newValue = value;
	} // end setValue

} // end public class DeviceBoxValueChangedEvent
