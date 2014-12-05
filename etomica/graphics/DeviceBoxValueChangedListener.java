/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

package etomica.graphics;

/**
 * 
 * Interface that should be implemented to receive a value changed event
 * when a DeviceBox updateValue method is called.
 *
 */
public interface DeviceBoxValueChangedListener extends java.util.EventListener {
	public void deviceBoxValueChanged(DeviceBoxValueChangedEvent ev);
} // end public interface DeviceBoxValueChangedListener
