import java.io.*;
import java.util.*;

/**
 * Main class for the Bayesian Network inference assignment.
 * Handles input parsing, network loading, and output formatting.
 */
public class Ex1 {
    /**
     * Main method that reads input, processes queries, and writes output.
     * 
     * @param args command line arguments (not used)
     * @throws Exception if an error occurs during execution
     */
    public static void main(String[] args) throws Exception {
        // Read input file
        BufferedReader reader = new BufferedReader(new FileReader("input_for_alarm2.txt"));
        String networkFile = reader.readLine();
        
        // Load the Bayesian network
        BayesianNetwork network = new BayesianNetwork();
        network.loadFromXML(networkFile);
        
        // Process queries and collect results
        List<String> results = new ArrayList<>();
        String line;
        while ((line = reader.readLine()) != null) {
            if (line.trim().isEmpty()) {
                continue;
            }
            
            BayesianNetwork.Result result = processQuery(network, line);
            results.add(result.toString());
        }
        reader.close();
        
        // Write results to output file
        BufferedWriter writer = new BufferedWriter(new FileWriter("output.txt"));
        for (int i = 0; i < results.size(); i++) {
            writer.write(results.get(i));
            if (i < results.size() - 1) {
                writer.newLine();
            }
        }
        writer.close();
    }
    
    /**
     * Processes a query and returns the result.
     * 
     * @param network the Bayesian network
     * @param query the query string from the input file
     * @return the result of the query
     */
    private static BayesianNetwork.Result processQuery(BayesianNetwork network, String query) {
        // Check if it's a joint probability query
        if (!query.contains("|")) {
            return processJointProbabilityQuery(network, query);
        } else {
            return processConditionalProbabilityQuery(network, query);
        }
    }
    
    /**
     * Processes a joint probability query.
     * 
     * @param network the Bayesian network
     * @param query the joint probability query string
     * @return the result of the query
     */
    private static BayesianNetwork.Result processJointProbabilityQuery(BayesianNetwork network, String query) {
        // Parse the query: P(B=F,E=T,A=T,M=T,J=F)
        String assignmentsStr = query.substring(query.indexOf('(') + 1, query.indexOf(')'));
        String[] assignmentParts = assignmentsStr.split(",");
        
        Map<String, String> assignments = new HashMap<>();
        for (String assignment : assignmentParts) {
            String[] parts = assignment.split("=");
            assignments.put(parts[0], parts[1]);
        }
        
        return network.calculateJointProbability(assignments);
    }
    
    /**
     * Processes a conditional probability query.
     * 
     * @param network the Bayesian network
     * @param query the conditional probability query string
     * @return the result of the query
     */
    private static BayesianNetwork.Result processConditionalProbabilityQuery(BayesianNetwork network, String query) {
        // Parse the query: P(B=T|J=T,M=T),1
        String[] queryParts = query.split(",(?=[0-9])");  // Split on comma followed by a number
        int algorithmNumber = Integer.parseInt(queryParts[1]);
        
        String probabilityQuery = queryParts[0];
        String[] conditionParts = probabilityQuery.split("\\|");
        
        // Parse query variable
        String queryAssignment = conditionParts[0].substring(conditionParts[0].indexOf('(') + 1);
        String[] queryParts2 = queryAssignment.split("=");
        Map.Entry<String, String> queryVar = new AbstractMap.SimpleEntry<>(queryParts2[0], queryParts2[1]);
        
        // Parse evidence variables
        Map<String, String> evidence = new HashMap<>();
        if (conditionParts.length > 1) {
            String evidenceStr = conditionParts[1].substring(0, conditionParts[1].indexOf(')'));
            String[] evidenceParts = evidenceStr.split(",");
            
            for (String evidencePart : evidenceParts) {
                String[] parts = evidencePart.split("=");
                evidence.put(parts[0], parts[1]);
            }
        }
        
        // Execute the appropriate algorithm
        switch (algorithmNumber) {
            case 1:
                return network.inferenceByEnumeration(queryVar, evidence);
            case 2:
                return network.variableElimination(queryVar, evidence);
            case 3:
                return network.algorithm3(queryVar, evidence);
            default:
                throw new IllegalArgumentException("Unknown algorithm number: " + algorithmNumber);
        }
    }
} 