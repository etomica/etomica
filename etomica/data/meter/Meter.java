package etomica.data.meter;

import etomica.EtomicaElement;
import etomica.data.DataSource;
import etomica.phase.PhaseDependent;


/**
 * A Phase-dependent DataSource.  Subclasses must implement the
 * getData(Phase) method
 *
 * @author David Kofke
 */

public interface Meter extends DataSource, EtomicaElement, PhaseDependent {

}
