import java.util.*;

// Represents a single node (random variable) in the Bayesian Network
class BayesianNode {
    private String name;                      // Node name
    private List<String> outcomes;            // Possible outcomes (values)
    private List<String> parents;             // Parent node names
    private Map<List<String>, Double> cpt;    // Conditional Probability Table (CPT)

    // Constructor initializes the node with a name and possible outcomes
    public BayesianNode(String name, List<String> outcomes) {
        this.name = name;
        this.outcomes = outcomes;
        this.parents = new ArrayList<>();
        this.cpt = new HashMap<>();
    }

    // Getters and setters
    public String getName() {
        return name;
    }

    public List<String> getOutcomes() {
        return outcomes;
    }

    public List<String> getParents() {
        return parents;
    }

    public void setParents(List<String> parents) {
        this.parents = parents;
    }

    public Map<List<String>, Double> getCPT() {
        return cpt;
    }

    public void setCPT(Map<List<String>, Double> cpt) {
        this.cpt = cpt;
    }
}

// Represents the entire Bayesian Network
public class BayesianNetwork {
    private Map<String, BayesianNode> nodes; // Map of node names to BayesianNode objects

    // Constructor initializes the network with an empty node map
    public BayesianNetwork() {
        nodes = new HashMap<>();
    }

    // Adds a new node to the network
    public void addNode(BayesianNode node) {
        nodes.put(node.getName(), node);
    }

    // Retrieves a node by its name
    public BayesianNode getNode(String nodeName) {
        return nodes.get(nodeName);
    }

    // Retrieves the entire map of nodes
    public Map<String, BayesianNode> getNodes() {
        return nodes;
    }
}