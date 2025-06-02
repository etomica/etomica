package etomica.GasMOP;
import etomica.space3d.Vector3D;
import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;


public class SharedPoints {
    public Map<Vector3D, ArrayList<ArrayList<Integer>>> findSharedPositions(Map<Integer, Map<Integer, Vector3D>> positionMap) {
        Map<Vector3D, ArrayList<ArrayList<Integer>>> sharedMap = new HashMap<>();

        for (Map.Entry<Integer, Map<Integer, Vector3D>> outer : positionMap.entrySet()) {
            int i = outer.getKey();
            for (Map.Entry<Integer, Vector3D> inner : outer.getValue().entrySet()) {
                int j = inner.getKey();
                Vector3D vec = inner.getValue();

                // Prepare [i, j] as a list
                ArrayList<Integer> pair = new ArrayList<>();
                pair.add(i);
                pair.add(j);

                // Insert into map
                sharedMap.computeIfAbsent(vec, k -> new ArrayList<>()).add(pair);
            }
        }

        // Keep only vectors shared by more than one atom pair
       /* Map<Vector3D, ArrayList<ArrayList<Integer>>> result = new HashMap<>();
        for (Map.Entry<Vector3D, ArrayList<ArrayList<Integer>>> entry : sharedMap.entrySet()) {
            if (entry.getValue().size() > 1) {
                result.put(entry.getKey(), entry.getValue());
            }
        }*/
        Map<RoundedVector3D, ArrayList<ArrayList<Integer>>> tempMap = new HashMap<>();
        Map<RoundedVector3D, Vector3D> representativeKeyMap = new HashMap<>();
        for (Map.Entry<Vector3D, ArrayList<ArrayList<Integer>>> entry : sharedMap.entrySet()) {
            Vector3D originalKey = entry.getKey();
            RoundedVector3D roundedKey = new RoundedVector3D(originalKey);

            tempMap.putIfAbsent(roundedKey, new ArrayList<>());
            tempMap.get(roundedKey).addAll(entry.getValue());

            // Keep first original Vector3D that maps to this rounded key
            representativeKeyMap.putIfAbsent(roundedKey, originalKey);
        }

        // Rebuild final map with original Vector3D keys
        Map<Vector3D, ArrayList<ArrayList<Integer>>> result = new HashMap<>();
        for (Map.Entry<RoundedVector3D, ArrayList<ArrayList<Integer>>> entry : tempMap.entrySet()) {
            Vector3D representativeKey = representativeKeyMap.get(entry.getKey());
            result.put(representativeKey, entry.getValue());
        }

        return result;

    }

    class RoundedVector3D {
        double x, y, z;

        public RoundedVector3D(Vector3D v) {
            this.x = round(v.getX(0));
            this.y = round(v.getX(1));
            this.z = round(v.getX(2));
        }

        private double round(double value) {
            return new BigDecimal(value).setScale(3, RoundingMode.HALF_UP).doubleValue();
        }

        @Override
        public boolean equals(Object o) {
            if (!(o instanceof RoundedVector3D)) return false;
            RoundedVector3D other = (RoundedVector3D) o;
            return Double.compare(x, other.x) == 0 &&
                    Double.compare(y, other.y) == 0 &&
                    Double.compare(z, other.z) == 0;
        }

        @Override
        public int hashCode() {
            return Objects.hash(x, y, z);
        }

    }




    public static void main(String[] args) {
        Map<Integer, Map<Integer, Vector3D>> positionMap = new HashMap<>();

        Vector3D v1 = new Vector3D(1.0, 2.0, 3.0);
        Vector3D v2 = new Vector3D(4.0, 5.0, 6.0);

        positionMap.computeIfAbsent(2, k -> new HashMap<>()).put(2, v1);
        positionMap.computeIfAbsent(3, k -> new HashMap<>()).put(4, v1);
        positionMap.computeIfAbsent(5, k -> new HashMap<>()).put(6, v1); // Now 3 atoms at v1
        positionMap.computeIfAbsent(7, k -> new HashMap<>()).put(8, v2);
        positionMap.computeIfAbsent(9, k -> new HashMap<>()).put(10, v2);
        SharedPoints sharedPoints = new SharedPoints();

        Map<Vector3D, ArrayList<ArrayList<Integer>>> shared = sharedPoints.findSharedPositions(positionMap);

        for (Map.Entry<Vector3D, ArrayList<ArrayList<Integer>>> entry : shared.entrySet()) {
            System.out.println("Position " + entry.getKey() + " shared by: " + entry.getValue());
        }
    }

}
