## Description of Instances

Two sets of benchmark instances are used in this paper to test the proposed algorithms:

1. **Monaci et al. (2013):**[[1]](#ref1) Instances for SMKP.
2. **Ryu and Park (2021):** [[2]](#ref2)Additional instances for SMKP, generated following their method.

All instances share the same format, organized according to the following rules:

- **First line:**
   Contains three values:
  - The number of items (*n*),
  - The bin capacity,
  - The confidence level.
- **Next *n* lines:**
   Each line corresponds to one item and contains three parameters:
  - Profit of the item (used only in SMKP),
  - Average weight of the item,
  - Weight variance of the item.

## Reference

<a id="ref1">[1]</a> Monaci, M., Pferschy, U., & Serafini, P. (2013). Exact solution of the robust knapsack problem. *Computers & operations research*, *40*(11), 2625-2631. (https://doi.org/10.1016/j.cor.2013.05.005).  

<a id="ref2">[2]</a> Ryu, J., & Park, S. (2021). Robust solutions for stochastic and distributionally robust chance-constrained binary knapsack problems. *arXiv preprint arXiv:2105.11875*. ([ https://doi.org/10.48550/arXiv.2105.11875](https://doi.org/10.48550/arXiv.2105.11875)).  
