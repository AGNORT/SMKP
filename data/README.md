## Instances

We used two sets of benchmark instances to test our exact method:

1. **Monaci.zip**: derived from the original instances by Monaci et al. (2013) [[1]](#ref1) for the Robust Knapsack Problem.
2. **RyuPark.zip**: generated according to the instructions provided by Ryu and Park (2021) [[2]](#ref2).

All instances share the same following format:

- First line – Three values:
  - The number of items (*n*);
  - The bin capacity;
  - The confidence level.
- Next *n* lines – Each line corresponds to one item and contains its three parameters:
  - Profit;
  - Nominal weight;
  - Submodular weight.

## References

<a id="ref1">[1]</a> Monaci, M., Pferschy, U., & Serafini, P. (2013). Exact solution of the robust knapsack problem. *Computers & operations research*, *40*(11), 2625-2631. (https://doi.org/10.1016/j.cor.2013.05.005).  

<a id="ref2">[2]</a> Ryu, J., & Park, S. (2021). Robust solutions for stochastic and distributionally robust chance-constrained binary knapsack problems. *arXiv preprint arXiv:2105.11875*. ([ https://doi.org/10.48550/arXiv.2105.11875](https://doi.org/10.48550/arXiv.2105.11875)).  
