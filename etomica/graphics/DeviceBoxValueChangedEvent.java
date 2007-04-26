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
