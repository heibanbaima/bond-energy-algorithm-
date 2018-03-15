import java.util.*;
import java.io.*;

public class Fragmentation {
    public static void main(String[] args) {
        Fragmentation f = new Fragmentation();
        f.go();

    }

    private void go() {
        Fragmentation f = new Fragmentation();
        Queries q = f.readInput();
        bondEnergy(q);
    }

    // not so important: this method reads the AF and AQ input files
    private Queries readInput() {
        Queries queries = new Queries();

        String af = readFile("AF.txt");
        StringTokenizer st1 = new StringTokenizer(af, "\n");
        int cnt = 0;
        while (st1.hasMoreTokens()) {
            String line = st1.nextToken();
            StringTokenizer st2 = new StringTokenizer(line, " ");
            int freq = 0;
            cnt++;
            Query q = new Query("q" + cnt);
            while (st2.hasMoreTokens()) {
                String num = st2.nextToken();
                freq += Integer.parseInt(num);
            }
            q.freq = freq;
            queries.add(q);
        }
        System.out.println(queries.toString());
        String aq = readFile("AQ.txt");
        StringTokenizer st3 = new StringTokenizer(aq, "\n");
        int queryCnt = 0;
        while (st3.hasMoreTokens()) {
            String line = st3.nextToken();
            StringTokenizer st4 = new StringTokenizer(line, " ");
            Query q = queries.queries.elementAt(queryCnt);
            int aCnt = 1;
            while (st4.hasMoreTokens()) {
                String num = st4.nextToken();
                if (num.equals("1") && aCnt > 0) {
                    q.attributes.add("A" + (aCnt));
                }
                aCnt++;
            }
            queryCnt++;
        }
        for (Query query : queries.queries) {
            System.out.println(query);
        }
        System.out.println("");

        return queries;
    }

    // not important: just reads a file
    public static String readFile(String filename) {
        String page = "";
        try {
            BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
            String line;
            while ((line = br.readLine()) != null) {
                page += line + "\n";
            }
            br.close();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return page;
    }


    // this is the implementation of the algorithm
    public void bondEnergy(Queries queries) {
        // the calculation of the affinity matrix is done in getMatrix
        AffinityMatrix m = queries.getMatrix();
        System.out.println(m);
        System.out.println("");
        Vector<String> attr = new Vector<String>(queries.getAllAttributes());
        LinkedList<String> remaining = new LinkedList<String>(attr);
        LinkedList<String> result = new LinkedList<String>();


        // pick first two attributes
        Random rand = new Random();
        result.add(remaining.remove(rand.nextInt(remaining.size())));
        result.add(remaining.remove(rand.nextInt(remaining.size())));

        // place remaining attributes
        System.out.println("Place attributes:");
        for (String newAttr : remaining) {
            System.out.println("place " + newAttr);
            // test new attribute
            int maxCont = 0;
            int pos = 0;
            for (int i = 0; i <= result.size(); i++) {
                String one = null;
                String two = null;
                if (i > 0) one = result.get(i - 1);
                if (i < result.size()) two = result.get(i);
                int cont = m.cont(one, newAttr, two);
                if (cont > maxCont) {
                    maxCont = cont;
                    pos = i;
                }
                System.out.println("contribution at pos " + i + " = " + cont);
            }
            result.add(pos, newAttr);
            System.out.println("attribute " + newAttr + " is placed at pos " + pos + ": " + result);
            System.out.println("");
        }

        System.out.println("resulting order: " + result);
        System.out.println("");

        // find fragmentation cut
        System.out.println("find fragments:");
        TreeSet<String> frag1;
        TreeSet<String> frag2;
        LinkedList<TreeSet<String>> f1 = new LinkedList<TreeSet<String>>();
        LinkedList<TreeSet<String>> f2 = new LinkedList<TreeSet<String>>();
        int maxSq = 0;

        // we have to test (n chooses 2) possible fragementations
        for (int i = 0; i < attr.size(); i++) {
            for (int j = (i + 1); j < attr.size(); j++) {
                frag2 = new TreeSet<String>(get(i + 1, j + 1, result));
                frag1 = new TreeSet<String>(result);
                frag1.removeAll(frag2);

                int acc1 = queries.accesses(frag1);
                int acc2 = queries.accesses(frag2);
                int acc3 = queries.accesses(frag1, frag2);
                int sq = acc1 * acc2 - acc3 * acc3;
                if (i == 0 && j == 1) maxSq = sq;
                if (sq > maxSq) {
                    maxSq = sq;
                    f1.clear();
                    f2.clear();
                }

                if (sq >= maxSq) {
                    f1.add(frag1);
                    f2.add(frag2);
                }
                System.out.println("split at " + frag1 + " | " + frag2);
                System.out.println("accesses frag1 alone: " + acc1);
                System.out.println("accesses frag2 alone: " + acc2);
                System.out.println("accesses frag1 and frag2: " + acc3);
                System.out.println("split quality = " + sq);
                System.out.println("");
            }
        }
        System.out.println("optimal split(s) (sq = " + maxSq + "):");
        for (int i = 0; i < f1.size(); i++) {
            System.out.println(f1.get(i) + " | " + f2.get(i));
        }

    }

    // given a list, this method return everything between the positions "from" and "to"
    public LinkedList<String> get(int from, int to, LinkedList<String> list) {
        LinkedList<String> retList = new LinkedList<String>();
        for (int i = from; i < to; i++) {
            retList.add(list.get(i));
        }
        return retList;
    }


}

// this class represents a query, it has an id, a query frequency (the sum over all sites), and a set of attributes
class Query implements Comparable<Query> {
    String id;
    int freq;
    TreeSet<String> attributes = new TreeSet<String>();

    public Query(String id) {
        this.id = id;
    }


    // returns freq if "att" is contained in the query
    public int access(String att) {
        if (attributes.contains(att)) return freq;
        else return 0;
    }

    // returns freq if both attributes are contained in the query
    public int access(String a1, String a2) {
        TreeSet<String> atts = new TreeSet<String>();
        atts.add(a1);
        atts.add(a2);
        if (attributes.containsAll(atts)) return freq;
        else return 0;
    }


    // return freq if the fragment covers this query
    public int exAccesses(Set<String> frag) {
        if (frag.containsAll(attributes)) return freq;
        else return 0;
    }

    // returns freq if both fragments are accessed by this query
    public int exAccesses(Set<String> frag1, Set<String> frag2) {
        boolean accessF1 = false;
        boolean accessF2 = false;
        for (String att : attributes) {
            if (frag1.contains(att)) accessF1 = true;
            if (frag2.contains(att)) accessF2 = true;
        }
        if (accessF1 && accessF2) return freq;
        else return 0;
    }

    public String toString() {
        String s = id + ": ";
        for (String att : attributes) {
            s += att + " ";
        }
        s += freq;
        return s;
    }

    public int compareTo(Query query) {
        return id.compareTo(query.id);
    }
}

// this class holds all queries
class Queries {
    Vector<Query> queries = new Vector<Query>();

    public void add(Query q) {
        queries.add(q);
    }

    public int affinity(String att) {
        return affinity(att, att);
    }

    // affinity between two attributes
    public int affinity(String a1, String a2) {
        int sum = 0;
        for (Query query : queries) {
            sum += query.access(a1, a2);
        }
        return sum;
    }

    public TreeSet<String> getAllAttributes() {
        TreeSet<String> set = new TreeSet<String>();
        for (Query query : queries) {
            set.addAll(query.attributes);
        }
        return set;
    }

    public AffinityMatrix getMatrix() {
        TreeSet<String> atts = new TreeSet<String>();
        for (Query q : queries) {
            atts.addAll(q.attributes);
        }
        AffinityMatrix m = new AffinityMatrix();
        for (String a1 : atts) {
            for (String a2 : atts) {
                m.put(a1, a2, affinity(a1, a2));
            }
        }
        return m;
    }

    // returns the sum of freq of all queries covered by "fragment"
    public int accesses(Set<String> fragment) {
        int sum = 0;
        for (Query query : queries) {
            sum += query.exAccesses(fragment);
        }
        return sum;
    }

    // returns the sum of freq of all queries that need both fragments
    public int accesses(Set<String> frag1, Set<String> frag2) {
        int sum = 0;
        for (Query query : queries) {
            sum += query.exAccesses(frag1, frag2);
        }
        return sum;
    }

    public String toString() {
        String s = "";
        for (Query query : queries) {
            s += query + "\n";
        }
        return s;
    }
}


class AffinityMatrix {
    TreeMap<String, TreeMap<String, Integer>> matrix = new TreeMap<String, TreeMap<String, Integer>>();

    public int get(String rowId, String colId) {
        TreeMap<String, Integer> row = matrix.get(rowId);
        if (row == null) return 0;
        Integer val = row.get(colId);
        if (val == null) return 0;
        return val;
    }

    // symmetric put of a val
    public void putSym(String rowId, String colId, int val) {
        put(rowId, colId, val);
        put(colId, rowId, val);
    }

    // put a val
    public void put(String rowId, String colId, int val) {
        TreeMap<String, Integer> row = matrix.get(rowId);
        if (row == null) {
            row = new TreeMap<String, Integer>();
            matrix.put(rowId, row);
        }
        row.put(colId, new Integer(val));
    }

    // bond between to attributes
    public int bond(String attr1, String attr2) {
        if (attr1 == null || attr2 == null) return 0;
        if (attr1.equals(attr2)) return 0;
        TreeMap<String, Integer> row1 = matrix.get(attr1);
        TreeMap<String, Integer> row2 = matrix.get(attr2);
        Iterator<Map.Entry<String, Integer>> i1 = row1.entrySet().iterator();
        Iterator<Map.Entry<String, Integer>> i2 = row2.entrySet().iterator();
        int sum = 0;
        while (i1.hasNext()) {
            int val1 = i1.next().getValue();
            int val2 = i2.next().getValue();
            sum += val1 * val2;
        }
        return sum;
    }

    // contribution of putting "middle" between "before" and "after"
    public int cont(String before, String middle, String after) {
        return bond(before, middle) + bond(middle, after) - bond(before, after);
    }

    // just prints the matrix
    public String toString() {
        int maxColLen = 3;
        TreeSet<String> attributes = new TreeSet<String>();
        for (Map.Entry<String, TreeMap<String, Integer>> e : matrix.entrySet()) {
            if (e.getKey().length() > maxColLen) maxColLen = e.getKey().length();
            attributes.add(e.getKey());
            for (String s : e.getValue().keySet()) {
                if (s.length() > maxColLen) maxColLen = s.length();
                attributes.add(s);
            }
        }
        maxColLen++;
        String firstLine = spaces(maxColLen);
        String s = "";
        for (String a : attributes) {
            firstLine += spaces(maxColLen - a.length() - 1) + a + " ";
            String line = a + spaces(maxColLen - a.length());
            TreeMap<String, Integer> row = matrix.get(a);
            if (row != null) {
                for (String colId : attributes) {
                    Integer valI = row.get(colId);
                    if (valI == null) {
                        line += spaces(maxColLen);
                    } else {
                        String valS = valI.toString();
                        line += spaces(maxColLen - valS.length() - 1) + valS + " ";
                    }
                }
            }
            s += line + "\n";
        }
        return firstLine + "\n" + s;
    }

    private String spaces(int amount) {
        String s = "";
        for (int i = 0; i < amount; i++) {
            s += " ";
        }
        return s;
    }
}