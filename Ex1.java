import org.w3c.dom.*;
import javax.xml.parsers.*;
import java.io.File;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.PrintWriter;

enum EliminationMode {
    FIXED_ORDER,
    HEURISTIC
}

// Main class responsible for loading and displaying the Bayesian network
public class Ex1 {

    // Loads Bayesian Network from XML file located at the given path
    public static BayesianNetwork loadNetworkFromXML(String path) throws Exception {
        BayesianNetwork network = new BayesianNetwork();

        File file = new File(path);
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();
        Document doc = builder.parse(file);
        doc.getDocumentElement().normalize();

        // Parse and load variables (nodes) into the network
        NodeList variables = doc.getElementsByTagName("VARIABLE");
        for (int i = 0; i < variables.getLength(); i++) {
            Element varElement = (Element) variables.item(i);
            String name = varElement.getElementsByTagName("NAME").item(0).getTextContent();

            NodeList outcomesNode = varElement.getElementsByTagName("OUTCOME");
            List<String> outcomes = new ArrayList<>();
            for (int j = 0; j < outcomesNode.getLength(); j++) {
                outcomes.add(outcomesNode.item(j).getTextContent());
            }

            BayesianNode node = new BayesianNode(name, outcomes);
            network.addNode(node);
        }

        // Parse definitions (CPTs and parent relations) for each node
        NodeList definitions = doc.getElementsByTagName("DEFINITION");
        for (int i = 0; i < definitions.getLength(); i++) {
            Element defElement = (Element) definitions.item(i);
            String nodeName = defElement.getElementsByTagName("FOR").item(0).getTextContent();

            // Retrieve parent node names
            NodeList givenNodes = defElement.getElementsByTagName("GIVEN");
            List<String> parents = new ArrayList<>();
            for (int j = 0; j < givenNodes.getLength(); j++) {
                String parent = givenNodes.item(j).getTextContent().trim();
                parents.add(parent);
            }

            // Retrieve CPT values (will be detailed later)
            String[] cptValues = defElement.getElementsByTagName("TABLE").item(0).getTextContent().trim().split("\\s+");

            // Set parents and initialize CPT (placeholder, to be filled later)
            BayesianNode node = network.getNode(nodeName);
            Collections.reverse(parents);
            node.setParents(parents);
            parseCPT(network, node, cptValues);
        }

        return network;
    }
    
    // Given a BayesianNode, parses its CPT from the given XML element.
    
    private static void parseCPT(BayesianNetwork network, BayesianNode node, String[] cptValues) {
    List<String> parents = node.getParents();
    List<String> outcomes = node.getOutcomes();
    
    // Map to store the CPT
    Map<List<String>, Double> cpt = new LinkedHashMap<>();

    // Generate all possible parent outcome combinations
    List<List<String>> parentCombinations = generateParentCombinations(network, parents);

    int index = 0;
    // Fill the CPT entries with the corresponding probabilities
    for (List<String> parentOutcome : parentCombinations) {
        for (String outcome : outcomes) {
            List<String> key = new ArrayList<>(parentOutcome);
            key.add(outcome);
            System.out.println("Added to CPT for " + node.getName() + ": " + key + " = " + cptValues[index]);
            cpt.put(key, Double.parseDouble(cptValues[index]));
            index++;
        }
    }

    node.setCPT(cpt);
}

// Helper method to generate all combinations of parent outcomes recursively.
 
private static List<List<String>> generateParentCombinations(BayesianNetwork network, List<String> parents) {
    List<List<String>> combinations = new ArrayList<>();
    
    if (parents.isEmpty()) {
        combinations.add(new ArrayList<>());
        return combinations;
    }

    String firstParent = parents.get(0);
    List<String> restParents = parents.subList(1, parents.size());
    
    for (String outcome : network.getNode(firstParent).getOutcomes()) {
        for (List<String> combination : generateParentCombinations(network, restParents)) {
            List<String> newCombination = new ArrayList<>();
            newCombination.add(outcome);
            newCombination.addAll(combination);
            combinations.add(newCombination);
        }
    }
    
    return combinations;
}
/**
 * Loads queries from the specified input file.
 * @param filePath The path to the queries input file.
 * @return List of queries as strings.
 * @throws Exception if file reading fails.
 */
private static List<String> loadQueriesFromFile(String filePath) throws Exception {
    List<String> queries = new ArrayList<>();

    BufferedReader reader = new BufferedReader(new FileReader(filePath));

    // Skip first line (XML filename already loaded)
    reader.readLine();

    String line;
    while ((line = reader.readLine()) != null) {
        line = line.trim();
        if (!line.isEmpty()) {
            queries.add(line);
        }
    }

    reader.close();
    return queries;
}

/**
 * Calculates the joint probability for given evidence using simple inference.
 * @param network The Bayesian network.
 * @param evidence A map containing variable names and their observed values.
 * @return The joint probability.
 */
public static double simpleInference(BayesianNetwork network, Map<String, String> evidence) {
    double probability = 1.0;

    for (String variable : evidence.keySet()) {
        BayesianNode node = network.getNode(variable);
        List<String> parentValues = new ArrayList<>();

        // Retrieve the values of the parents from the evidence
        for (String parent : node.getParents()) {
            if (!evidence.containsKey(parent)) {
                throw new IllegalArgumentException("Missing parent value in evidence: " + parent);
            }
            parentValues.add(evidence.get(parent));
        }

        // Append the variable's own value to the key
        parentValues.add(evidence.get(variable));

        Double prob = node.getCPT().get(parentValues);
        if (prob == null) {
            throw new IllegalStateException("Missing CPT entry for: " + parentValues);
        }
        System.out.println("Looking up CPT for node: " + variable);
        System.out.println("Key used: " + parentValues);
        System.out.println("→ Found probability: " + prob);
        System.out.println("Current product: " + probability + " * " + prob + " = " + (probability * prob));
        System.out.println("----------------------");
        probability *= prob;
    }

    return probability;
}

private static String pickMinFactorVariable(Set<String> varsToCheck, List<List<String>> factorVars, BayesianNetwork net) {
    String bestVar = null;
    int minSize = Integer.MAX_VALUE;

    for (String var : varsToCheck) {
        Set<String> involved = new HashSet<>();
        for (List<String> vars : factorVars) {
            if (vars.contains(var)) involved.addAll(vars);
        }

        int size = 1;
        for (String v : involved) {
            size *= net.getNode(v).getOutcomes().size();
        }

        if (size < minSize) {
            minSize = size;
            bestVar = var;
        }
    }

    return bestVar;
}
public static double variableEliminationFixed(BayesianNetwork net, String queryVar, String queryVal, Map<String, String> evidence, List<String> order) {
    return runVariableElimination(net, queryVar, queryVal, evidence, EliminationMode.FIXED_ORDER, order);
}

public static double variableEliminationHeuristic(BayesianNetwork net, String queryVar, String queryVal, Map<String, String> evidence) {
    return runVariableElimination(net, queryVar, queryVal, evidence, EliminationMode.HEURISTIC, new ArrayList<>());
}

private static double runVariableElimination(
    BayesianNetwork network,
    String queryVar,
    String queryValue,
    Map<String, String> evidence,
    EliminationMode mode,
    List<String> fixedOrder
) {
    List<Map<List<String>, Double>> factors = new ArrayList<>();
    List<List<String>> factorVars = new ArrayList<>();

    // Step 1: Initialize factors from the CPTs of the network
    for (BayesianNode node : network.getNodes().values()) {
        List<String> vars = new ArrayList<>(node.getParents());
        vars.add(node.getName());
        factors.add(new LinkedHashMap<>(node.getCPT()));
        factorVars.add(vars);
    }

    // Step 2: Reduce each factor using the given evidence
    for (int i = 0; i < factors.size(); i++) {
        Map<List<String>, Double> factor = factors.get(i);
        List<String> vars = factorVars.get(i);
        Map<List<String>, Double> reduced = new LinkedHashMap<>();

        for (Map.Entry<List<String>, Double> entry : factor.entrySet()) {
            List<String> key = entry.getKey();
            boolean match = true;
            for (int j = 0; j < vars.size(); j++) {
                String var = vars.get(j);
                if (evidence.containsKey(var) && !evidence.get(var).equals(key.get(j))) {
                    match = false;
                    break;
                }
            }
            if (match) {
                reduced.put(key, entry.getValue());
            }
        }

        factors.set(i, reduced);
    }

    // Step 3: Identify variables to eliminate
    Set<String> toEliminate = new HashSet<>();
    for (String var : network.getNodes().keySet()) {
        if (!var.equals(queryVar) && !evidence.containsKey(var)) {
            toEliminate.add(var);
        }
    }

    int elimIndex = 0;
    while (!toEliminate.isEmpty()) {
        String elim;
        if (mode == EliminationMode.FIXED_ORDER) {
            elim = fixedOrder.get(elimIndex++);
        } else {
            elim = pickMinFactorVariable(toEliminate, factorVars, network);
        }
        toEliminate.remove(elim);

        // Collect all factors that mention the current variable
        List<Map<List<String>, Double>> toJoin = new ArrayList<>();
        List<List<String>> toJoinVars = new ArrayList<>();
        List<Integer> indexes = new ArrayList<>();

        for (int i = 0; i < factorVars.size(); i++) {
            if (factorVars.get(i).contains(elim)) {
                indexes.add(i);
                toJoin.add(factors.get(i));
                toJoinVars.add(factorVars.get(i));
            }
        }

        if (toJoin.isEmpty()) continue;

        // Remove these factors from the lists
        for (int i = indexes.size() - 1; i >= 0; i--) {
            int idx = indexes.get(i);
            factors.remove(idx);
            factorVars.remove(idx);
        }

        // Compute the product (join)
        List<String> allVars = new ArrayList<>();
        for (List<String> vars : toJoinVars) {
            for (String v : vars) {
                if (!allVars.contains(v)) {
                    allVars.add(v);
                }
            }
        }

        Map<List<String>, Double> joined = new LinkedHashMap<>();
        for (List<String> assignment : generateCombinations(network, allVars, evidence)) {
            double product = 1.0;
            for (int i = 0; i < toJoin.size(); i++) {
                List<String> vars = toJoinVars.get(i);
                List<String> key = new ArrayList<>();
                for (String v : vars) {
                    key.add(assignment.get(allVars.indexOf(v)));
                }

                Double val = toJoin.get(i).get(key);
                if (val == null) {
                    product = 0.0;
                    break;
                }
                product *= val;
            }
            joined.put(assignment, product);
        }

        // Sum out the eliminated variable
        Map<List<String>, Double> summed = new LinkedHashMap<>();
        List<String> remainingVars = new ArrayList<>(allVars);
        remainingVars.remove(elim);

        for (Map.Entry<List<String>, Double> entry : joined.entrySet()) {
            List<String> newKey = new ArrayList<>();
            for (String v : remainingVars) {
                newKey.add(entry.getKey().get(allVars.indexOf(v)));
            }

            summed.put(newKey, summed.getOrDefault(newKey, 0.0) + entry.getValue());
        }

        factors.add(summed);
        factorVars.add(remainingVars);
    }

    // Step 4: Multiply remaining factors
    Map<List<String>, Double> finalF = factors.get(0);
    List<String> finalVars = factorVars.get(0);

    for (int i = 1; i < factors.size(); i++) {
        Map<List<String>, Double> f2 = factors.get(i);
        List<String> vars2 = factorVars.get(i);

        Set<String> allSet = new LinkedHashSet<>(finalVars);
        allSet.addAll(vars2);
        List<String> allVars = new ArrayList<>(allSet);

        Map<List<String>, Double> newF = new LinkedHashMap<>();
        for (List<String> assignment : generateCombinations(network, allVars, evidence)) {
            List<String> key1 = new ArrayList<>();
            List<String> key2 = new ArrayList<>();
            for (String v : finalVars) key1.add(assignment.get(allVars.indexOf(v)));
            for (String v : vars2) key2.add(assignment.get(allVars.indexOf(v)));

            Double v1 = finalF.get(key1);
            Double v2 = f2.get(key2);
            if (v1 != null && v2 != null) {
                newF.put(assignment, v1 * v2);
            }
        }

        finalF = newF;
        finalVars = allVars;
    }

    // Step 5: Normalize the distribution
    double total = 0.0, desired = 0.0;
    int idx = finalVars.indexOf(queryVar);

    for (Map.Entry<List<String>, Double> entry : finalF.entrySet()) {
        total += entry.getValue();
        if (entry.getKey().get(idx).equals(queryValue)) {
            desired += entry.getValue();
        }
    }

    return desired / total;
}
private static List<List<String>> generateCombinations(BayesianNetwork network, List<String> vars, Map<String, String> evidence) {
    List<List<String>> res = new ArrayList<>();
    backtrack(res, new ArrayList<>(), vars, 0, network, evidence);
    return res;
}

private static void backtrack(List<List<String>> res, List<String> curr, List<String> vars, int i, BayesianNetwork network, Map<String, String> evidence) {
    if (i == vars.size()) {
        res.add(new ArrayList<>(curr));
        return;
    }
    String var = vars.get(i);
    if (evidence.containsKey(var)) {
        curr.add(evidence.get(var));
        backtrack(res, curr, vars, i + 1, network, evidence);
        curr.remove(curr.size() - 1);
    } else {
        for (String val : network.getNode(var).getOutcomes()) {
            curr.add(val);
            backtrack(res, curr, vars, i + 1, network, evidence);
            curr.remove(curr.size() - 1);
        }
    }
}

public static double bruteForceConditional(BayesianNetwork network, String queryVar, String queryValue, Map<String, String> evidence) {
    List<String> hiddenVars = new ArrayList<>();
    for (String var : network.getNodes().keySet()) {
        if (!evidence.containsKey(var) && !var.equals(queryVar)) {
            hiddenVars.add(var);
        }
    }

    double numerator = 0.0;
    double denominator = 0.0;

    // numerator: queryVar = queryValue
    Map<String, String> withQuery = new HashMap<>(evidence);
    withQuery.put(queryVar, queryValue);

    List<Map<String, String>> allNumeratorAssignments = generateAssignments(network, hiddenVars);
    for (Map<String, String> assign : allNumeratorAssignments) {
        Map<String, String> full = new HashMap<>(withQuery);
        full.putAll(assign);
        numerator += simpleInference(network, full);
    }

    // denominator: sum over all possible values of queryVar
    for (String val : network.getNode(queryVar).getOutcomes()) {
        Map<String, String> ev = new HashMap<>(evidence);
        ev.put(queryVar, val);
        List<Map<String, String>> assignments = generateAssignments(network, hiddenVars);
        for (Map<String, String> assign : assignments) {
            Map<String, String> full = new HashMap<>(ev);
            full.putAll(assign);
            denominator += simpleInference(network, full);
        }
    }

    return numerator / denominator;
}

private static List<Map<String, String>> generateAssignments(BayesianNetwork network, List<String> vars) {
    List<Map<String, String>> result = new ArrayList<>();
    backtrackAssign(result, new HashMap<>(), vars, 0, network);
    return result;
}

private static void backtrackAssign(List<Map<String, String>> result, Map<String, String> current, List<String> vars, int index, BayesianNetwork network) {
    if (index == vars.size()) {
        result.add(new HashMap<>(current));
        return;
    }
    String var = vars.get(index);
    for (String val : network.getNode(var).getOutcomes()) {
        current.put(var, val);
        backtrackAssign(result, current, vars, index + 1, network);
        current.remove(var);
    }
}

public static void main(String[] args) {
    try {
        // Load the Bayesian Network
        BayesianNetwork network = loadNetworkFromXML("alarm_net.xml");
        System.out.println("Bayesian Network loaded successfully!");
        System.out.println("Nodes loaded into the network:");

        // Display structure: names, outcomes, parents
        for (String nodeName : network.getNodes().keySet()) {
            BayesianNode node = network.getNode(nodeName);
            System.out.println("Node: " + node.getName());
            System.out.println("Outcomes: " + node.getOutcomes());
            System.out.println("Parents: " + node.getParents());
            System.out.println("-----------------------------------");
        }

        // Display each node's CPT for manual inspection
        System.out.println("\n--- CPT Verification ---");
        for (BayesianNode node : network.getNodes().values()) {
            System.out.println("Node: " + node.getName());
            System.out.println("CPT entries:");
            for (Map.Entry<List<String>, Double> entry : node.getCPT().entrySet()) {
                System.out.println(entry.getKey() + " = " + entry.getValue());
            }
            System.out.println("------------------------");
        }

        // Load queries file (will be used later)
        List<String> queries = loadQueriesFromFile("input.txt");
        System.out.println("\nQueries loaded:");
        for (String query : queries) {
            System.out.println(query);
        }

        // ----------------------------------------
        // Example evidence and queries
        // ----------------------------------------
        Map<String, String> fullEvidence = new HashMap<>();
        fullEvidence.put("B", "F");
        fullEvidence.put("E", "T");
        fullEvidence.put("A", "T");
        fullEvidence.put("M", "T");
        fullEvidence.put("J", "F");

        Map<String, String> partialEvidence = new HashMap<>();
        partialEvidence.put("J", "T");
        partialEvidence.put("M", "T");

        String queryVar = "B";
        String queryVal = "T";
        List<String> fixedOrder = Arrays.asList("E", "A");

        // ----------------------------------------
        // Run and compare all inference methods
        // ----------------------------------------

        System.out.println("\n--- Inference Results ---");

        double simple = simpleInference(network, fullEvidence);
        System.out.println("Simple Inference:\nP(B=F, E=T, A=T, M=T, J=F) = " + simple);

        double fixed = variableEliminationFixed(network, queryVar, queryVal, partialEvidence, fixedOrder);
        System.out.println("\nVariable Elimination (Fixed Order):\nP(B=T | J=T, M=T) = " + fixed);

        double heuristic = variableEliminationHeuristic(network, queryVar, queryVal, partialEvidence);
        System.out.println("\nVariable Elimination (Heuristic Order):\nP(B=T | J=T, M=T) = " + heuristic);

        double brute = bruteForceConditional(network, queryVar, queryVal, partialEvidence);
        System.out.println("\nBrute Force Conditional:\nP(B=T | J=T, M=T) = " + brute);

        // ----------------------------------------
        // Sanity checks / cross-verification
        // ----------------------------------------

        System.out.println("\n--- Verification ---");
        if (Math.abs(fixed - brute) < 1e-6) {
            System.out.println("✅ Fixed Order VE matches Brute Force");
        } else {
            System.out.println("❌ Fixed Order VE does NOT match Brute Force");
        }

        if (Math.abs(heuristic - brute) < 1e-6) {
            System.out.println("✅ Heuristic VE matches Brute Force");
        } else {
            System.out.println("❌ Heuristic VE does NOT match Brute Force");
        }

        // Optional: Output results to file
        PrintWriter writer = new PrintWriter("output.txt", "UTF-8");
        writer.println("P(B=F, E=T, A=T, M=T, J=F) = " + simple);
        writer.println("P(B=T | J=T, M=T) [Fixed Order] = " + fixed);
        writer.println("P(B=T | J=T, M=T) [Heuristic Order] = " + heuristic);
        writer.println("P(B=T | J=T, M=T) [Brute Force] = " + brute);
        writer.close();

    } catch (Exception e) {
        System.err.println("❌ Failed to load or execute inference: " + e.getMessage());
        e.printStackTrace();
    }
}
}