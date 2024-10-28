# Multi Index Sets

```@example misets
using SparseGridsKit

n, k = 2, 1
mi_set = create_smolyak_miset(n, k)
```

```@example misets
mi_set_new = MISet([[1,3]]) 
combined_mi_set = add_mi(mi_set, mi_set_new)
```

```@example misets
margin_set = get_margin(mi_set)
```

```@example misets
reduced_margin_set = get_reduced_margin(mi_set)
```