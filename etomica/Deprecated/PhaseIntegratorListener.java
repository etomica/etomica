package simulate;
import java.util.*;

public interface PhaseIntegratorListener extends EventListener {
    public void phaseIntegratorNotify(PhaseIntegratorEvent evt);
}