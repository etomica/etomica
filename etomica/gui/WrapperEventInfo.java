package etomica.gui;

import java.lang.reflect.Method;
import java.beans.EventSetDescriptor;

/**
 * Event listener hookup information.
 *
 * This has the same information as WrapperEventTarget but in a more
 * accessible way.  It is quite arguable that we should just extend that
 * class but the final insight came very late in the Nov BDK cycle.
 * Will revisit in the future
 */

public class WrapperEventInfo {

    Object targetBean;
    String adaptorClassName;
    String eventSetName;

    /**
     * Represents an adaptor connection between a source (implicitly given) and
     * a single action mediated by the adaptor on a target.
     * 
     * There are two types of adaptors, those for property bound changes, and
     * more generic ones for actions.
     *
     * @param target - the target object
     * @param klass - the adaptor class represented by this hookup
     * @param esn - the EventSetName of this source
     */
    public WrapperEventInfo(Object target, String hookupClassName,
			    String esn) {
	targetBean = target;
	adaptorClassName = hookupClassName;
	eventSetName = esn;
    }

    // The target for the hookup
    public Object getTargetBean() {
	return targetBean;
    }

    // The adaptor Class
    public String getAdaptorClassName() {
	return adaptorClassName;
    }

    // The target Interface type
    public String getEventSetName() {
	return eventSetName;
    }
}
