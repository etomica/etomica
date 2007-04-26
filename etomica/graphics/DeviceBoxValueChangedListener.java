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
