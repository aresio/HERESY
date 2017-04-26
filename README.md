# HERESY

HERESY (Highly Efficient REaction SYstem simulator) is a tool for the simplified simulation of Reaction Systems (RS). HERESY provides the user with an intuitive Graphic User Interface and it is able to leverage the GPU to offload reactions' calculations.

# Using HERESY

To run HERESY use the following command:

`python heresy.py`

and then use the GUI to add a reaction system, along with context to be investigated, formalized as follows.

_Reactions:_ reactions are entered as single lines containing the sets of reactants, inhibitors, and products. The elements of the sets are separated by spaces; the three sets are separated by commas. E.g., a reaction ({r1, r2}, {i1}, {p1, p2, p3}) is entered as:

`r1 r2, i1, p1 p2 p3`

_Context:_ the set of chemicals for each iteration is entered as a single line separated by spaces. E.g., a context {s1, s2, s3} is entered as:

`s1 s2 s3`

# Citing HERESY

TBA

