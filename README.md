# Bayesian Network Inference Engine

## 🔍 Project Overview

This project implements an inference engine for Bayesian networks.  
It supports three primary algorithms for answering probabilistic queries:

- **Simple Inference** (Joint probability computation)
- **Variable Elimination** with fixed or heuristic variable ordering
- **Brute Force Conditional Probability**

All queries are parsed from an input file and evaluated against a Bayesian network defined in XML format.

---

## 🧠 Features

- Parses `.xml` Bayesian network structure (nodes, outcomes, CPTs)
- Supports queries in the format:  
  `P(A=T,B=F)` or `P(A=T|B=F,C=T),<algorithm_number>`
- Outputs results to `output.txt` in the format:  
  ```
  <probability>,<additions>,<multiplications>
  ```
  with 5 decimal places of precision
- Tracks computational effort (additions and multiplications)
- Includes internal debugging utilities to print CPTs and network structure

---

## ⚙️ Algorithms

### 1. Simple Inference (`method = 1`)
Computes the **joint probability** of a complete assignment by traversing the CPTs of all nodes and multiplying the relevant probabilities.

### 2. Variable Elimination (Fixed / Heuristic)
Eliminates hidden variables using factor multiplication and marginalization.

- **Fixed Order:** Uses a predetermined list of variables to eliminate.
- **Heuristic Order:** Selects the next variable to eliminate based on minimal intermediate factor size (min-fill heuristic).

### 3. Brute Force Conditional
Calculates:
```
P(Query | Evidence) = P(Query ∧ Evidence) / P(Evidence)
```
by enumerating all hidden variable assignments and summing over them.

---

## 🧪 Input File Structure

The input file (`input.txt`) has the following format:

```
<path_to_network.xml>
<query1>
<query2>
...
```

Examples:
```
alarm_net.xml
P(B=F,E=T,A=T,M=T,J=F)
P(B=T|J=T,M=T),1
P(J=T|B=T),2
```

---

## 📝 Output Format

The output file (`output.txt`) contains one line per query:
```
0.07375,12,38
```
- **0.07375** – The resulting probability  
- **12** – Number of additions (during normalization)  
- **38** – Number of multiplications used in computation

---

## 🏗️ Project Structure

- `Ex1.java` – Main entry point, XML parser, inference algorithms
- `BayesianNetwork.java` – Data model for nodes and the network
- `input.txt` – Contains network filename and queries
- `alarm_net.xml` – Sample Bayesian network definition
- `output.txt` – Output results (generated after execution)

---

## 🚀 How to Run

1. Compile:
   ```bash
   javac Ex1.java
   ```

2. Run:
   ```bash
   java Ex1
   ```

3. View the output:
   ```
   cat output.txt
   ```

---

## 📚 Example Network: `alarm_net.xml`

Includes variables:
- `E, B, A, J, M`
with CPTs defined for:
- Burglary, Earthquake, Alarm, JohnCalls, MaryCalls

---

## 👨‍💻 Authors
- Developed by Yoav Peretz as part of the **AI Algorithms course**  
  B.Sc. in Computer Science

---

## 📦 Future Improvements
- Add graphical visualization of the network
- Integrate external Bayesian network formats (e.g., BIFXML)
- Optimize VE with caching and pruning