# GraphRandomisation

Control over condition ordering is achieved using Graph Theory, as outlined in: *Brooks JL (2012). Counterbalancing for serial order carryover effects in experimental conditions orders. Psychological Methods, 17(4): 600-614*. In fact, the Python code I’ve written is adapted from the MATLAB scripts he provided. So please cite that.

## Using the module

To use the module in your code, you first need to import it – which additionally means the file needs to be placed either somewhere in the Python path or in the same folder as your experiment script.

```
import GraphRandomisation
```

Or if you don’t want to write ```GraphRandomisation.``` in front of any function you call from it:
```
from GraphRandomisation import *
```
The next step is to create a dictionary of sequence information. The sequence information minimally needs to have a key ‘conditions’ that specifies either the number of conditions or is the full list of every condition. You can also specify the ‘order’ to control for (default is 1); the number of ‘repetitions’ of each combination (default is 1); whether to start from a particular node (default is None, which results in a random choice); and whether to ‘omit_self_adjacencies’, which is a Boolean of whether to prevent each condition from occurring next to itself (default is False). Note that omitting self-adjacencies should only ever be used with 1st-degree order. Using it with higher degrees will produce weird results. In addition, the start variable corresponds to which node to start from, rather than which condition. This means that for 2nd-degree order and higher, this specifies a particular tuple of conditions to start (and end) with. For example, with 2nd-degree ordering the 1st node is (0,0), the 2nd node is (0,1), and so on.

Here’s an example of a sequence specification of 5 conditions, control over 1st-degree order, with 5 repetitions, and starting with the 1st condition (N.B. Python starts counting from 0). Note that when controlling for 1st-degree order each condition will occur a number of times equal to the number of conditions per repetition. So in this example, each condition would occur 20 times (5 (conditions) x 4 (repetitions)):
```
sequence_info = {'conditions':5, 'order':1, 'repetitions':4, 'start':0}
```
The sequence can then be used to create a trial sequence:
```
sequence = GraphRandomisation.generate_sequence(sequence_info)
```
If you want to try and produce a more uniform trial sequence you can use the following function:
```
sequence = GraphRandomisation.generate_uniform_sequence(sequence_info, n_sequences=50)
```
This will generate a number of sequences equal to the number specified (if omitted, it will default to 100). It will then return the sequence that is the most uniform from amongst them (using the procedure outlined in Sohn et al., 1997).
The returned sequence will simply be a list of numbers from 0...k-1 (where k is the number of conditions – remember Python starts counting from 0). Here’s an example of a sequence produced using 4 conditions (with an order of 1, and only 1 repetition):
```
[3, 3, 0, 1, 1, 0, 3, 1, 2, 1, 3, 2, 0, 0, 2, 2, 3]
```
Note that the sequence is 17 items in length: each of the 4 conditions occurs 4 times, plus an extra trial of condition 3 at the beginning of the sequence. This extra trial should be dropped from the analysis.
The formula for determining the length of the final sequence is:

$`r * k^{(c+1)} + c`$

Where r is the number of repeats, k is the number of conditions, and c is the degree of order controlled for. The extra c trials are added to the beginning of the sequence as padding, and should be dropped from any analysis. This means the full trial list can only be multiples of: $`k^{(c+1)}`$ (plus the padding: c). If omitting self-adjacencies, the formula changes slightly to:

$`r * k * (k-1)^c + c`$

There are functions available to provide the number of trials per repetition for some specified sequence information (as detailed above), as well as the total number of trials that would be produced in the sequence:
```
trials_per_repeat = GraphRandomisation.get_trials_per_repeat(sequence_info)
total_trials_in_sequence = GraphRandomisation.get_trials_in_sequence(sequence_info)
```
If you’d like to use unequal trial orders, you can specify your own Adjacency Matrix. An adjacency matrix is a 2D matrix that specifies the edges of the graph, with the exits from nodes as columns and the entries into nodes as rows (or vice versa, I can’t quite remember). For example, an entry in the first row and second column would specify an edge travelling from the 1st node to the 2nd node. The adjacency matrices I create from the sequence information you give me are uniform – every entry in the matrix has the same value. However, if you know what you’re doing you could create an imbalanced matrix (remember rows and columns relate to nodes, not conditions) and pass that into the following function (also passing in a list of the nodes that correspond to each row/column):
```
sequence = GraphRandomisation.generate_sequence_using_matrix(matrix,nodes)
```
The one restriction is that in order to be able to take a Eulerian Circuit through the graph, there needs to be an equal number of entry and exits into each node. Practically, this means that the sum of the ith column must be equal to the sum of the ith row (i.e, it needs to be symmetrical around the diagonal).


